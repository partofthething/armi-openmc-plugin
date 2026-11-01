"""
Microscopic multigroup cross sections from OpenMC, written as ISOTXS.

The model
---------
A representative block is modeled alone as an infinite lattice (its footprint with reflective
sides) or, when its cross section group names a driver, at the center of a colorset of driver
blocks:

* Cartesian: a 3x3 lattice with the target in the middle and reflective sides, which is an exact
  infinite checkerboard of the colorset.
* Hexagonal: the target with one ring of driver blocks, inside a cylinder of homogenized driver
  with a white boundary. Hexagonal colorsets do not tile, so this is an approximation of an
  infinite driver region around one target.

A block with no fissionable material and no driver is run as a fixed-source problem with a
fission spectrum source, which gives it a slowing-down spectrum.

Reduction
---------
For every nuclide ``i`` the tallies give reaction rates ``R_i,g`` over the target block, and the
flux tally gives the volume-integrated group flux ``Phi_g``. The flux-volume weighted microscopic
cross section, the one that reproduces the block's reaction rates when multiplied by the
block-averaged number density, is ``sigma_i,g = R_i,g / (N_i Phi_g)``.

* Reaction cross sections (total, fission, capture, (n,2n), (n,a), (n,p), (n,d), (n,t)) and
  nu, the fission energy and the group velocities come from collision (expected-value)
  estimators, which score every nuclide at every collision and so converge for every nuclide in
  the block, trace ones included, without the cost track-length scoring has in a pin lattice.
* The scattering matrix takes its magnitude from the expected scattering rate and its shape
  (where neutrons go, and how many come out) from an analog estimator. The matrix is
  production-weighted, (n,2n) and friends included, and goes into ISOTXS as the one scattering
  block (elastic, P0). ARMI's absorption counts (n,2n) as a loss, so the neutron balance closes.
  A nuclide with too few analog collisions to shape its matrix (a trace nuclide) gets its
  scattering in-group.
* Transport cross sections subtract the analog P1 scattering moment from the total.
* Fission spectra come from analog outgoing-energy tallies; a nuclide with no analog fissions
  gets the block's fission spectrum.

Nuclides that a block lacks but ARMI wants in the library are added to every material at a trace
density (``openmcXSTraceDensity``), which leaves the transport unchanged while letting the
collision estimators see them.
"""

import json
import math
import os
from dataclasses import dataclass, field

import h5py
import numpy as np
import openmc
from armi import runLog
from armi.nucDirectory import nuclideBases
from armi.nuclearDataIO import xsCollections, xsLibraries, xsNuclides
from armi.nuclearDataIO.cccc import isotxs
from armi.reactor import grids
from armi.reactor.blocks import HexBlock
from scipy import sparse

from armi_openmc import blocks as blockModels
from armi_openmc import shapes
from armi_openmc.core import makeSettings
from armi_openmc.energy import energyFilterBounds, toArmiOrder, upperBounds
from armi_openmc.materials import KELVIN_OFFSET, MaterialFactory, NuclideMapper
from armi_openmc.options import OpenMCOptions

EV_TO_J = 1.602176634e-19
REACTIONS = [
    "total",
    "fission",
    "nu-fission",
    "kappa-fission",
    "scatter",
    "(n,gamma)",
    "(n,2n)",
    "(n,a)",
    "(n,p)",
    "(n,d)",
    "(n,t)",
]
ISOTXS_EXTRAS = {"(n,a)": "nalph", "(n,p)": "np", "(n,2n)": "n2n", "(n,d)": "nd", "(n,t)": "nt"}
ELASTIC_P0 = 100
#: Estimator of the flux and reaction-rate tallies. The collision (expected-value) estimator scores
#: every nuclide at every collision, so trace nuclides converge with the flux, but costs nothing on
#: the many short track segments of a pin lattice, where a track-length tally of dozens of
#: nuclides and reactions would do a cross section lookup per nuclide per segment.
REACTION_ESTIMATOR = "collision"


@dataclass
class XSModel:
    """An OpenMC model for one cross section group, with what is needed to reduce its tallies."""

    xsID: str
    block: object
    model: openmc.Model
    #: ARMI nuclide -> block-averaged number density used in the reduction (trace included).
    densities: dict
    #: ARMI nuclide -> [(OpenMC nuclide, fraction)].
    expansion: dict
    omcNuclides: list
    energyBounds: np.ndarray
    temperatureK: float
    eigenvalue: bool = True
    notes: list = field(default_factory=list)


class XSModelBuilder:
    """Builds the OpenMC model of one representative block (and its driver)."""

    def __init__(self, options: OpenMCOptions, hexGrid=None):
        self.options = options
        self.mapper = NuclideMapper(options.crossSectionsPath)
        self.hexOrientation = blockModels.hexOrientationOf(hexGrid) if hexGrid is not None else "x"
        self.hexGrid = hexGrid

    def build(self, block, xsID, nuclides, driver=None) -> XSModel:
        """
        Parameters
        ----------
        block : Block
            The representative block.
        xsID : str
            Two-letter cross section ID the nuclides are labeled with.
        nuclides : iterable of str
            ARMI nuclide names the library must contain for this group.
        driver : Block, optional
            Block to surround the target with.
        """
        opts = self.options
        expansion = {n: self.mapper.expand(n) for n in sorted(set(nuclides))}
        expansion = {n: e for n, e in expansion.items() if e}
        trace = {n: opts.xsTraceDensity for n in expansion}
        materials = MaterialFactory(self.mapper, opts.thermalScattering, traceNuclides=trace)
        builder = blockModels.BlockBuilder(materials, self.hexOrientation)
        heterogeneous = opts.blockModel == "heterogeneous"
        target = builder.build(block, heterogeneous)

        isHex = isinstance(block, HexBlock)
        outline = blockModels.blockOutline(block, self.hexOrientation)
        notes = []
        if driver is not None:
            driverModel = builder.build(driver, heterogeneous)
            if isHex:
                root, targetCell, bounds = self._hexColorset(
                    block, target, driverModel, outline, builder
                )
            else:
                root, targetCell, bounds = self._rectColorset(target, driverModel, outline)
            notes.append(f"colorset with driver {driver.getMicroSuffix()}")
            fissionable = bool(target.fissionableCells or driverModel.fissionableCells)
        else:
            region = outline.interior(boundaryType="reflective")
            targetCell = openmc.Cell(name=f"{xsID} target", fill=target.universe, region=region)
            root = [targetCell]
            r = outline.circumscribedRadius
            bounds = (-r, -r, r, r)
            fissionable = bool(target.fissionableCells)
            notes.append("infinite lattice")

        bottom = openmc.ZPlane(-0.5, boundary_type="reflective")
        top = openmc.ZPlane(0.5, boundary_type="reflective")
        for cell in root:
            cell.region = (
                (cell.region & +bottom & -top) if cell.region is not None else (+bottom & -top)
            )
        geometry_ = openmc.Geometry(openmc.Universe(cells=root))
        mats = openmc.Materials(materials.materials)
        mats.cross_sections = opts.crossSectionsPath

        x0, y0, x1, y1 = bounds
        box = (x0, y0, -0.5, x1, y1, 0.5)
        if fissionable:
            settings = makeSettings(opts, box)
        else:
            notes.append("fixed source (no fissionable material)")
            source = openmc.IndependentSource(
                space=openmc.stats.Box(box[:3], box[3:]),
                energy=openmc.stats.Watt(a=0.988e6, b=2.249e-6),
            )
            settings = makeSettings(opts, box, runMode="fixed source", source=source)
            settings.batches = max(1, opts.batches - opts.inactive)

        densities = {}
        for n in expansion:
            densities[n] = block.getNumberDensity(n) or opts.xsTraceDensity
        omcNuclides = sorted({name for e in expansion.values() for name, _ in e})
        bounds_ = energyFilterBounds(opts.groupStructure)
        tallies = _xsTallies(targetCell, omcNuclides, bounds_)
        model = openmc.Model(geometry_, mats, settings, tallies)
        return XSModel(
            xsID=xsID,
            block=block,
            model=model,
            densities=densities,
            expansion=expansion,
            omcNuclides=omcNuclides,
            energyBounds=bounds_,
            temperatureK=block.getAverageTempInC() + KELVIN_OFFSET,
            eigenvalue=fissionable,
            notes=notes,
        )

    def _rectColorset(self, target, driver, outline):
        px, py = outline.size
        holder = openmc.Cell(name="target", fill=target.universe)
        targetUniverse = openmc.Universe(cells=[holder])
        lattice = openmc.RectLattice(name="colorset")
        lattice.pitch = (px, py)
        lattice.lower_left = (-1.5 * px, -1.5 * py)
        d = driver.universe
        lattice.universes = [[d, d, d], [d, targetUniverse, d], [d, d, d]]
        box = shapes.Outline("rectangle", (3 * px, 3 * py))
        root = openmc.Cell(
            name="colorset", fill=lattice, region=box.interior(boundaryType="reflective")
        )
        return [root], holder, (-1.5 * px, -1.5 * py, 1.5 * px, 1.5 * py)

    def _hexColorset(self, block, target, driver, outline, builder):
        pitch = outline.size[0]
        holder = openmc.Cell(name="target", fill=target.universe)
        targetUniverse = openmc.Universe(cells=[holder])
        driverBlock = driver.block
        smeared = builder.materials.make(
            f"{driverBlock.getName()} smeared",
            driverBlock.getNumberDensities(),
            driverBlock.getAverageTempInC(),
        )
        outer = openmc.Universe(cells=[openmc.Cell(name="driver region", fill=smeared)])
        lattice = openmc.HexLattice(name="colorset")
        lattice.orientation = "y" if self.hexOrientation == "x" else "x"
        lattice.center = (0.0, 0.0)
        lattice.pitch = (pitch,)
        lattice.universes = [[driver.universe] * 6, [targetUniverse]]
        lattice.outer = outer
        # A cylinder holding the area of three rings (19 blocks): two rings modeled in detail,
        # the third smeared.
        radius = math.sqrt(19 * outline.area / math.pi)
        region = -openmc.ZCylinder(r=radius, boundary_type="white")
        root = openmc.Cell(name="colorset", fill=lattice, region=region)
        return [root], holder, (-radius, -radius, radius, radius)


def _xsTallies(targetCell, nuclides, bounds):
    domain = openmc.CellFilter([targetCell])
    energy = openmc.EnergyFilter(bounds)
    energyout = openmc.EnergyoutFilter(bounds)

    flux = openmc.Tally(name="xs flux")
    flux.filters = [domain, energy]
    flux.scores = ["flux", "inverse-velocity"]
    flux.estimator = REACTION_ESTIMATOR

    reactions = openmc.Tally(name="xs reactions")
    reactions.filters = [domain, energy]
    reactions.nuclides = nuclides
    reactions.scores = REACTIONS
    reactions.estimator = REACTION_ESTIMATOR

    matrix = openmc.Tally(name="xs scatter matrix")
    matrix.filters = [domain, energy, energyout]
    matrix.nuclides = nuclides
    matrix.scores = ["scatter", "nu-scatter"]
    matrix.estimator = "analog"

    p1 = openmc.Tally(name="xs scatter p1")
    p1.filters = [domain, energy, openmc.LegendreFilter(1)]
    p1.nuclides = nuclides
    p1.scores = ["scatter"]
    p1.estimator = "analog"

    chi = openmc.Tally(name="xs chi")
    chi.filters = [domain, energyout]
    chi.nuclides = nuclides + ["total"]
    chi.scores = ["nu-fission"]
    chi.estimator = "analog"
    return openmc.Tallies([flux, reactions, matrix, p1, chi])


# ---------------------------------------------------------------------------------------------
@dataclass
class GroupConstants:
    """Microscopic group constants of one ARMI nuclide, in ARMI (fast-first) group order."""

    total: np.ndarray
    transport: np.ndarray
    fission: np.ndarray
    nu: np.ndarray
    chi: np.ndarray
    capture: np.ndarray
    extras: dict  # ISOTXS name -> array
    scatter: np.ndarray  # [to, from], production-weighted
    efiss: float  # J per fission
    density: float


def reduce(xsModel: XSModel, statepointPath):
    """
    Microscopic group constants for every nuclide of an XS model, from its statepoint.

    Returns
    -------
    (dict, dict)
        {ARMI nuclide: GroupConstants} and run information (keff, group velocities, flux).
    """
    with openmc.StatePoint(statepointPath) as sp:
        keff = sp.keff if xsModel.eigenvalue else None
        fluxT = sp.get_tally(name="xs flux")
        ng = len(xsModel.energyBounds) - 1
        flux = fluxT.get_values(scores=["flux"]).ravel()
        invVel = fluxT.get_values(scores=["inverse-velocity"]).ravel()
        rx = sp.get_tally(name="xs reactions")
        mat = sp.get_tally(name="xs scatter matrix")
        p1 = sp.get_tally(name="xs scatter p1")
        chiT = sp.get_tally(name="xs chi")

        def rxValues(nuc, score):
            return rx.get_values(scores=[score], nuclides=[nuc]).ravel()

        perOmc = {}
        omcDens = _openmcDensities(xsModel)
        for nuc in xsModel.omcNuclides:
            dens = omcDens[nuc]
            denom = np.where(flux > 0.0, dens * flux, np.inf)
            values = {score: rxValues(nuc, score) / denom for score in REACTIONS}
            fis = rxValues(nuc, "fission")
            values["nu"] = np.divide(
                rxValues(nuc, "nu-fission"), fis, out=np.zeros(ng), where=fis > 0
            )
            values["efissEV"] = (
                rxValues(nuc, "kappa-fission").sum() / fis.sum() if fis.sum() > 0 else 0.0
            )
            scat = mat.get_values(scores=["scatter"], nuclides=[nuc]).reshape(ng, ng)
            nuScat = mat.get_values(scores=["nu-scatter"], nuclides=[nuc]).reshape(ng, ng)
            values["matrix"] = _scatterMatrix(values["scatter"], scat, nuScat, values["(n,2n)"])
            moments = p1.get_values(scores=["scatter"], nuclides=[nuc]).reshape(ng, 2)
            values["p1"] = moments[:, 1] / denom
            values["chi"] = chiT.get_values(scores=["nu-fission"], nuclides=[nuc]).ravel()
            perOmc[nuc] = values
        blockChi = chiT.get_values(scores=["nu-fission"], nuclides=["total"]).ravel()

    blockChi = blockChi / blockChi.sum() if blockChi.sum() > 0 else np.zeros(ng)
    empty = toArmiOrder(flux <= 0.0)
    if empty.any():
        runLog.info(
            f"OpenMC XS {xsModel.xsID}: {int(empty.sum())} of {ng} groups saw no flux; they get "
            "the cross sections of the nearest group that did, scattering in-group"
        )
    constants = {}
    for armiName, parts in xsModel.expansion.items():
        gc = _combine(parts, perOmc, blockChi, xsModel.densities[armiName], ng)
        constants[armiName] = _fillEmptyGroups(gc, empty)
    velocity = np.divide(flux, invVel, out=np.zeros(ng), where=invVel > 0)
    velocity = np.where(velocity > 0, velocity, _midpointSpeeds(xsModel.energyBounds))
    info = {
        "keff": keff,
        "velocity": toArmiOrder(velocity),
        "flux": toArmiOrder(flux),
    }
    return constants, info


def _nearestPopulated(empty):
    """For every group, the index of the nearest group that is not empty."""
    populated = np.flatnonzero(~empty)
    if populated.size == 0:
        return np.arange(len(empty))
    return np.array([populated[np.argmin(np.abs(populated - g))] for g in range(len(empty))])


def _fillEmptyGroups(gc, empty):
    """
    Give groups that saw no flux the constants of the nearest group that did.

    A fast spectrum leaves the thermal groups of a broad structure empty, and a library of zeros
    there would make those groups singular in any solver that carries them. The neutrons that
    could reach them are too few to matter, so any finite, physical value will do; the nearest
    group's, with its scattering kept in-group, is the least surprising one.
    """
    if not empty.any():
        return gc
    source = _nearestPopulated(empty)
    for name in ("total", "transport", "fission", "nu", "capture"):
        values = getattr(gc, name)
        values[empty] = values[source[empty]]
    for values in gc.extras.values():
        values[empty] = values[source[empty]]
    outScatter = gc.scatter.sum(axis=0)  # [to, from]: total scattering out of each source group
    for g in np.flatnonzero(empty):
        gc.scatter[:, g] = 0.0
        gc.scatter[g, g] = outScatter[source[g]]
    return gc


def _midpointSpeeds(bounds):
    """Neutron speed (cm/s) at the logarithmic midpoint of each group, in ascending energy."""
    mid = np.sqrt(bounds[:-1] * bounds[1:])  # eV
    massEV = 939.56542e6  # neutron rest energy, eV
    return 2.99792458e10 * np.sqrt(2.0 * mid / massEV)


def _openmcDensities(xsModel):
    """Block-averaged density of each OpenMC nuclide, built from the ARMI ones it comes from."""
    out = {}
    for armiName, parts in xsModel.expansion.items():
        for name, frac in parts:
            out[name] = out.get(name, 0.0) + xsModel.densities[armiName] * frac
    return out


def _scatterMatrix(scatterTL, scatterAnalog, nuScatterAnalog, n2nTL):
    """
    Production-weighted scattering matrix [from, to] per unit flux and density.

    Magnitude from the expected (collision) scattering rate, shape and multiplicity from the analog
    tallies; in-group when the analog tallies saw nothing.
    """
    ng = len(scatterTL)
    matrix = np.zeros((ng, ng))
    for g in range(ng):
        events = scatterAnalog[g].sum()
        produced = nuScatterAnalog[g].sum()
        if events > 0.0 and produced > 0.0:
            matrix[g] = scatterTL[g] * (produced / events) * nuScatterAnalog[g] / produced
        else:
            matrix[g, g] = scatterTL[g] + n2nTL[g]
    return matrix


def _combine(parts, perOmc, blockChi, density, ng):
    """Sum an ARMI nuclide's constituent OpenMC nuclides (an element's isotopes) by abundance."""

    def weighted(key):
        return sum(frac * perOmc[name][key] for name, frac in parts)

    fission = weighted("fission")
    nuFission = weighted("nu-fission")
    nu = np.divide(nuFission, fission, out=np.zeros(ng), where=fission > 0)
    chi = sum(perOmc[name]["chi"] for name, _ in parts)
    if fission.sum() > 0:
        chi = chi / chi.sum() if chi.sum() > 0 else blockChi.copy()
    else:
        chi = np.zeros(ng)
    efissWeights = [frac * perOmc[name]["fission"].sum() for name, frac in parts]
    efiss = (
        sum(w * perOmc[name]["efissEV"] for w, (name, _) in zip(efissWeights, parts, strict=True))
        / sum(efissWeights)
        if sum(efissWeights) > 0
        else 0.0
    )
    total = weighted("total")
    transport = total - weighted("p1")
    matrix = weighted("matrix")  # [from, to], OpenMC ascending energy order
    return GroupConstants(
        total=toArmiOrder(total),
        transport=toArmiOrder(transport),
        fission=toArmiOrder(fission),
        nu=toArmiOrder(nu),
        chi=toArmiOrder(chi),
        capture=toArmiOrder(weighted("(n,gamma)")),
        extras={key: toArmiOrder(weighted(score)) for score, key in ISOTXS_EXTRAS.items()},
        # flip both axes into fast-first order, then transpose to ARMI's [to, from]
        scatter=np.flip(matrix, axis=(0, 1)).T.copy(),
        efiss=efiss * EV_TO_J,
        density=density,
    )


# ---------------------------------------------------------------------------------------------
def captureEnergy(mapper: NuclideMapper, armiName):
    """Q value of radiative capture in J, from the HDF5 data (abundance-weighted for elements)."""
    total = 0.0
    for name, frac in mapper.expand(armiName):
        path = mapper.dataPaths.get(name)
        if not path:
            continue
        try:
            with h5py.File(path, "r") as f:
                group = f[name]["reactions"].get("reaction_102")
                q = float(group.attrs["Q_value"]) if group is not None else 0.0
        except (OSError, KeyError):
            q = 0.0
        total += frac * q
    return total * EV_TO_J


def buildLibrary(results, groupStructure, mapper: NuclideMapper):
    """
    An ARMI :py:class:`~armi.nuclearDataIO.xsLibraries.IsotxsLibrary` from reduced constants.

    Parameters
    ----------
    results : list of (XSModel, dict, dict)
        Each model with its constants and run information, from :py:func:`reduce`.
    """
    lib = xsLibraries.IsotxsLibrary()
    upper = upperBounds(groupStructure)
    ng = len(upper)
    lib.neutronEnergyUpperBounds = upper
    velocities = [info["velocity"] for _, _, info in results if np.any(info["velocity"] > 0)]
    lib.neutronVelocity = np.mean(velocities, axis=0) if velocities else np.zeros(ng)

    #: The group flux each XS ID was weighted with, for comparisons (not part of ISOTXS).
    lib.openmcFlux = {xsModel.xsID: info["flux"] for xsModel, _, info in results}
    maxUp = maxDown = 0
    for xsModel, constants, _info in results:
        for armiName, gc in constants.items():
            nuc, up, down = _xsNuclide(lib, xsModel, armiName, gc, mapper, ng)
            maxUp, maxDown = max(maxUp, up), max(maxDown, down)
            lib[nuc.containerKey] = nuc

    meta = lib.isotxsMetadata
    meta["label"] = "ISOTXS"
    meta["fileId"] = 1
    meta["numGroups"] = ng
    meta["maxUpScatterGroups"] = maxUp
    meta["maxDownScatterGroups"] = maxDown
    meta["maxScatteringOrder"] = 0
    meta["fileWideChiFlag"] = 0
    meta["maxScatteringBlocks"] = 1
    meta["subblockingControl"] = 1
    meta["libraryLabel"] = "OpenMC"
    meta["chi"] = None
    meta["minimumNeutronEnergy"] = 1.0e-5
    return lib


def _xsNuclide(lib, xsModel, armiName, gc: GroupConstants, mapper, ng):
    base = nuclideBases.byName[armiName]
    label = base.label + xsModel.xsID
    nuc = xsNuclides.XSNuclide(lib, label)
    nuc._base = base
    micros = nuc.micros
    micros.numGroups = ng
    micros.total = gc.total.reshape(ng, 1)
    micros.transport = gc.transport.reshape(ng, 1)
    micros.nGamma = gc.capture
    fissile = bool(np.any(gc.fission > 0))
    micros.fission = gc.fission
    micros.neutronsPerFission = gc.nu
    micros.chi = gc.chi if fissile else xsCollections.XSCollection.getDefaultXs(ng)
    for key, values in gc.extras.items():
        setattr(micros, key, values)
    micros.strpd = xsCollections.XSCollection.getDefaultXs(ng)

    scatter = gc.scatter
    jband, jj, up, down = _bandStructure(scatter)
    micros.elasticScatter = sparse.csr_matrix(scatter)
    micros.inelasticScatter = sparse.csr_matrix((ng, ng))
    micros.n2nScatter = sparse.csr_matrix((ng, ng))

    md = nuc.isotxsMetadata
    md["nuclideId"] = armiName
    md["libName"] = "OPENMC"
    md["isoIdent"] = base.label
    md["amass"] = float(base.weight)
    md["efiss"] = float(gc.efiss)
    md["ecapt"] = float(captureEnergy(mapper, armiName))
    md["temp"] = float(xsModel.temperatureK)
    md["sigPot"] = 0.0
    md["adens"] = float(gc.density)
    md["classif"] = 0
    md["chiFlag"] = 1 if fissile else 0
    md["fisFlag"] = 1 if fissile else 0
    for key in ISOTXS_EXTRAS.values():
        md[key] = 1
    md["ltot"] = 1
    md["ltrn"] = 1
    md["strpd"] = 0
    md["scatFlag"] = np.array([ELASTIC_P0])
    md["ords"] = np.array([1])
    md["jband"] = {(g, 0): jband[g] for g in range(ng)}
    md["jj"] = {(g, 0): jj[g] for g in range(ng)}
    return nuc, up, down


def _bandStructure(scatter):
    """
    ISOTXS band storage of a [to, from] matrix: for each destination group, the number of source
    groups stored (JBAND) and the position of the in-group term within them (JJ).
    """
    ng = scatter.shape[0]
    jband, jj = [], []
    up = down = 0
    for g in range(ng):
        cols = np.nonzero(scatter[g])[0]
        lo = min(cols.min(), g) if cols.size else g
        hi = max(cols.max(), g) if cols.size else g
        jband.append(int(hi - lo + 1))
        jj.append(int(hi - g + 1))
        up = max(up, hi - g)
        down = max(down, g - lo)
    return jband, jj, int(up), int(down)


def writeLibrary(lib, path):
    """Write the ISOTXS, and beside it ``<path>.flux.json``: each XS ID's OpenMC group flux."""
    isotxs.writeBinary(lib, path)
    flux = getattr(lib, "openmcFlux", None)
    if flux:
        with open(f"{path}.flux.json", "w") as f:
            json.dump({k: [float(v) for v in values] for k, values in flux.items()}, f, indent=1)
    runLog.info(f"Wrote {len(lib)} nuclides for XS IDs {sorted(lib.xsIDs)} to {path}")


def blockHexGrid(block):
    """The hex grid a block sits in, for the orientation of its hexagons."""
    try:
        grid = block.parent.spatialLocator.grid
    except AttributeError:
        return None
    return grid if isinstance(grid, grids.HexGrid) else None


def runDirFor(options, xsID, suffix=""):
    return os.path.join(options.runDir, f"xs-{xsID}{suffix}")
