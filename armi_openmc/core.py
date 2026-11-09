"""
A whole-core OpenMC model of an ARMI reactor.

The core's spatial grid becomes an OpenMC lattice (hexagonal or rectangular) of assembly
universes; each assembly universe stacks its block universes between the block elevations ARMI
holds. Each assembly gets its own universes and materials, so every ARMI block can be tallied on
its own and later changed on its own.

Symmetry
--------
* A full core is modeled as is.
* A quarter Cartesian core with reflective boundaries is modeled as the quarter it is, cut by
  reflective planes through the origin (the ``native`` treatment).
* Any other partial core (third-core hex, periodic or eighth Cartesian) is unfolded to full core:
  each assembly's universe is placed again at its symmetric equivalents. Tallies on a block then
  sum over its copies, and the mapping back divides by the number of copies.
"""

import math
from dataclasses import dataclass, field

import numpy as np
import openmc
from armi import runLog
from armi.reactor import geometry, grids

from armi_openmc import blocks as blockModels
from armi_openmc.energy import energyFilterBounds
from armi_openmc.materials import MaterialFactory, NuclideMapper
from armi_openmc.multigroup import MultigroupMaterialFactory
from armi_openmc.options import OpenMCOptions

ROUND = 9
#: The multigroup library of a multigroup-mode model, beside its model.xml.
MGXS_FILE = "mgxs.h5"


def multigroupLibrary(core, options: OpenMCOptions):
    """The ISOTXS of the multigroup mode: ``openmcMultigroupLibrary``, else the core's own."""
    if options.multigroupLibrary:
        from armi.nuclearDataIO.cccc import isotxs

        return isotxs.readBinary(options.multigroupLibrary)
    if core.lib is None:
        raise ValueError(
            "OpenMC multigroup mode needs a cross section library: run a cross section kernel "
            "(genXS with an xsKernel) or set openmcMultigroupLibrary to an ISOTXS"
        )
    return core.lib


@dataclass
class CoreModel:
    """An OpenMC model of a core plus what is needed to map results back onto it."""

    model: openmc.Model
    #: Cell filled by each ARMI block's universe, keyed by the block.
    blockCells: dict = field(default_factory=dict)
    #: How many times each block appears in the model (symmetric copies of its assembly).
    copies: dict = field(default_factory=dict)
    #: True when the model is a reflective cut of the core (ARMI's partial blocks as they are).
    cut: bool = False
    #: Energy bounds of the multigroup tallies, ascending in eV.
    energyBounds: np.ndarray | None = None
    blockModels: dict = field(default_factory=dict)
    #: Fraction of the full core the modeled domain represents (1 for full or unfolded cores).
    domainFraction: float = 1.0
    #: The multigroup data of a multigroup-mode model, which the executer writes as MGXS_FILE.
    multigroup: object = None


class CoreModelBuilder:
    """Builds a :py:class:`CoreModel` from an ARMI core."""

    def __init__(self, core, options: OpenMCOptions):
        self.core = core
        self.options = options
        self.grid = core.spatialGrid
        self.isHex = isinstance(self.grid, grids.HexGrid)
        if options.multigroup:
            self.mapper = None
            self.materials = MultigroupMaterialFactory(
                multigroupLibrary(core, options),
                transportCorrected=options.multigroupScattering == "transport-corrected",
            )
        else:
            self.mapper = NuclideMapper(options.crossSectionsPath)
            self.materials = MaterialFactory(self.mapper, options.thermalScattering)
        self.blocks = blockModels.BlockBuilder(
            self.materials,
            blockModels.hexOrientationOf(self.grid) if self.isHex else "x",
            options=options,
        )
        self._zPlanes = {}
        self._void = openmc.Universe(name="void", cells=[openmc.Cell(name="void")])

    # ----------------------------------------------------------------------------------------
    def build(self) -> CoreModel:
        opts = self.options
        spherical = [b for b in self.core.iterBlocks() if blockModels.isSpherical(b)]
        if spherical:
            raise ValueError(
                f"{spherical[0]} is a spherical unit cell, which has no place in a core lattice. "
                "Run it with openmc-xs, or model the pebble bed with prismatic blocks."
            )
        result = CoreModel(model=None)
        positions, result.cut = self._positions()
        result.domainFraction = 1.0 / self.core.powerMultiplier if result.cut else 1.0

        assemblyUniverses = {}
        for a in self.core:
            assemblyUniverses[a] = self._assemblyUniverse(a, result)
        for a, places in positions.items():
            for b in result.blockCells:
                if b.parent is a:
                    result.copies[b] = len(places)

        lattice, extent = self._coreLattice(positions, assemblyUniverses)
        region, bounds = self._boundingRegion(extent, result.cut)
        root = openmc.Cell(name="core", fill=lattice, region=region)
        geometry_ = openmc.Geometry(openmc.Universe(name="root", cells=[root]))

        materials = openmc.Materials(self.materials.materials)
        settings = self._settings(bounds, self._fissionableBounds(result, positions, bounds))
        if opts.multigroup:
            materials.cross_sections = (
                MGXS_FILE  # the executer writes it beside model.xml, with its full path
            )
            settings.energy_mode = "multi-group"
            result.energyBounds = self.materials.bounds
            result.multigroup = self.materials
            self.materials.report()
        else:
            materials.cross_sections = opts.crossSectionsPath
            result.energyBounds = energyFilterBounds(opts.groupStructure)
        if self.mapper is not None:
            self.mapper.reportMissing()
        tallies = self._tallies(result)
        result.model = openmc.Model(geometry_, materials, settings, tallies)
        if opts.plotGeometry:
            result.model.plots = self._plots(bounds)
        nHet = sum(1 for m in result.blockModels.values() if m.heterogeneous)
        runLog.info(
            f"OpenMC core model: {len(result.blockCells)} blocks ({nHet} heterogeneous) in "
            f"{len(positions)} assemblies, {len(materials)} materials, "
            f"{'reflective cut' if result.cut else 'full core'}"
        )
        return result

    # ----------------------------------------------------------------------------------------
    def _positions(self):
        """Lattice indices (i, j) each assembly occupies, and whether the model is a cut."""
        sym = self.core.symmetry
        positions = {a: [tuple(a.spatialLocator.getCompleteIndices()[:2])] for a in self.core}
        if sym.domain == geometry.DomainType.FULL_CORE:
            return positions, False
        if (
            self.options.symmetry == "native"
            and not self.isHex
            and sym.domain == geometry.DomainType.QUARTER_CORE
            and sym.boundary == geometry.BoundaryType.REFLECTIVE
        ):
            return positions, True
        runLog.info(f"OpenMC: unfolding the {sym} core to full core")
        for places in positions.values():
            i, j = places[0]
            for eq in self.grid.getSymmetricEquivalents((i, j)):
                eq = tuple(eq[:2])
                if eq not in places:
                    places.append(eq)
        return positions, False

    def _heterogeneous(self, b):
        if self.options.multigroup or self.options.blockModel != "heterogeneous":
            return False
        flags = self.options.heterogeneousFlags
        if not flags:
            return True
        from armi.reactor.flags import Flags

        return any(b.hasFlags(Flags.fromStringIgnoreErrors(f)) for f in flags)

    def _zPlane(self, z):
        key = round(float(z), ROUND)
        if key not in self._zPlanes:
            self._zPlanes[key] = openmc.ZPlane(z0=key)
        return self._zPlanes[key]

    def _assemblyUniverse(self, a, result):
        cells = []
        sliceZ = self.options.axialSlice
        for b in a:
            if sliceZ is not None and not (b.p.zbottom <= sliceZ < b.p.ztop):
                continue
            model = self.blocks.build(b, self._heterogeneous(b), zOrigin=b.p.zbottom)
            result.blockModels[b] = model
            if sliceZ is not None:
                region = None
            else:
                region = +self._zPlane(b.p.zbottom) & -self._zPlane(b.p.ztop)
            cell = openmc.Cell(name=f"{b.getName()}", fill=model.universe, region=region)
            cells.append(cell)
            result.blockCells[b] = cell
        if sliceZ is None and len(a):
            cells.append(
                openmc.Cell(name=f"{a.getName()} below", region=-self._zPlane(a[0].p.zbottom))
            )
            cells.append(
                openmc.Cell(name=f"{a.getName()} above", region=+self._zPlane(a[-1].p.ztop))
            )
        if not cells:
            return self._void
        return openmc.Universe(name=a.getName(), cells=cells)

    # ----------------------------------------------------------------------------------------
    def _coreLattice(self, positions, universes):
        # Lattice positions without an assembly hold whatever lies outside the assemblies.
        self._empty = self._outerUniverse()
        coords = {}
        for a, places in positions.items():
            for ij in places:
                x, y, _ = self.grid[ij[0], ij[1], 0].getLocalCoordinates()
                coords[(float(x), float(y))] = universes[a]
        if self.isHex:
            lattice = self._hexLattice(positions)
            half = self.grid.pitch / math.sqrt(3.0)  # circumradius of an assembly
            radius = max(math.hypot(x, y) for x, y in coords) + half
            if self.options.coreRadius:
                if self.options.coreRadius < radius:
                    raise ValueError(
                        f"openmcCoreRadius {self.options.coreRadius} cm cuts through the "
                        f"assemblies, which reach {radius:.2f} cm"
                    )
                radius = self.options.coreRadius
            extent = ("circle", radius)
        else:
            px, py = (float(p) for p in self.grid.pitch)
            xs = [c[0] for c in coords]
            ys = [c[1] for c in coords]
            lattice = openmc.RectLattice(name="core")
            lattice.pitch = (px, py)
            lattice.lower_left = (min(xs) - px / 2.0, min(ys) - py / 2.0)
            nx = int(round((max(xs) - min(xs)) / px)) + 1
            ny = int(round((max(ys) - min(ys)) / py)) + 1
            lattice.universes = np.full((ny, nx), self._empty, dtype=object)
            extent = ("box", min(xs) - px / 2, max(xs) + px / 2, min(ys) - py / 2, max(ys) + py / 2)
            if self.options.coreRadius:
                # a cylinder around the square lattice, e.g. assemblies in a pool or a drum
                reach = max(math.hypot(abs(x) + px / 2, abs(y) + py / 2) for x, y in coords)
                if self.options.coreRadius < reach:
                    raise ValueError(
                        f"openmcCoreRadius {self.options.coreRadius} cm cuts through the "
                        f"assemblies, whose corners reach {reach:.2f} cm"
                    )
                extent = ("circle", self.options.coreRadius)
        lattice.outer = self._empty
        grid = lattice.universes
        for (x, y), universe in coords.items():
            idx = lattice.find_element((x, y, 0.0))[0]
            blockModels._setNested(grid, lattice.get_universe_index(idx), universe)
        lattice.universes = grid
        return lattice, extent

    def _outerUniverse(self):
        """What lies outside the assembly lattice: void, or the ``outerFill`` block's mixture."""
        if not self.options.outerFill:
            return self._void
        block = next(
            (b for b in self.core.iterBlocks() if b.getType() == self.options.outerFill), None
        )
        if block is None:
            # an assembly design that is not loaded in the core, such as the water of a pool
            blueprints = getattr(self.core.parent, "blueprints", None)
            designs = blueprints.assemblies.values() if blueprints is not None else []
            block = next(
                (b for a in designs for b in a if b.getType() == self.options.outerFill), None
            )
        if block is None:
            raise ValueError(
                f"openmcOuterFill: no block of type {self.options.outerFill} in the core "
                "or in the blueprints' assembly designs"
            )
        mat = self.materials.make(
            f"outside the assemblies ({block.getType()})",
            block.getNumberDensities(),
            block.getAverageTempInC(),
        )
        return openmc.Universe(name="outer fill", cells=[openmc.Cell(name="outer fill", fill=mat)])

    def _hexLattice(self, positions):
        rings = max(self.grid.getRingPos(ij)[0] for places in positions.values() for ij in places)
        x0, y0, _ = self.grid[0, 0, 0].getLocalCoordinates()
        lattice = openmc.HexLattice(name="core")
        lattice.orientation = blockModels.hexLatticeOrientationOf(self.grid)
        lattice.center = (float(x0), float(y0))
        lattice.pitch = (float(self.grid.pitch),)
        lattice.universes = [[self._empty] * max(1, 6 * (r - 1)) for r in range(rings, 0, -1)]
        return lattice

    def _boundingRegion(self, extent, cut):
        opts = self.options
        radial = opts.radialBoundary
        if extent[0] == "circle":
            if radial == "reflective":
                raise ValueError(
                    "A hexagonal core cannot have a reflective radial boundary in OpenMC."
                )
            radius = extent[1]
            region = -openmc.ZCylinder(r=radius, boundary_type=radial)
            xy = (-radius, radius, -radius, radius)
        else:
            _, x0, x1, y0, y1 = extent
            if cut:
                x0, y0 = max(x0, 0.0), max(y0, 0.0)
            lo = "reflective" if cut else radial
            region = (
                +openmc.XPlane(x0, boundary_type=lo)
                & -openmc.XPlane(x1, boundary_type=radial)
                & +openmc.YPlane(y0, boundary_type=lo)
                & -openmc.YPlane(y1, boundary_type=radial)
            )
            xy = (x0, x1, y0, y1)
        if opts.axialSlice is not None:
            z0, z1 = opts.axialSlice - 0.5, opts.axialSlice + 0.5
            bottom = openmc.ZPlane(z0, boundary_type="reflective")
            top = openmc.ZPlane(z1, boundary_type="reflective")
        else:
            z0 = min(a[0].p.zbottom for a in self.core)
            z1 = max(a[-1].p.ztop for a in self.core)
            bottom, top = self._zPlane(z0), self._zPlane(z1)
            bottom.boundary_type = opts.bottomBoundary
            top.boundary_type = opts.topBoundary
        region &= +bottom & -top
        return region, (xy[0], xy[2], z0, xy[1], xy[3], z1)

    # ----------------------------------------------------------------------------------------
    def _settings(self, bounds, sourceBounds):
        settings = makeSettings(self.options, bounds, sourceBounds=sourceBounds)
        return settings

    def _fissionableBounds(self, result, positions, bounds):
        """The box around every block holding fissionable material, clipped to the model."""
        half = self.grid.pitch / 2.0 if self.isHex else max(self.grid.pitch) / 2.0
        if self.isHex:
            half *= 2.0 / math.sqrt(3.0)
        boxes = []
        for b, model in result.blockModels.items():
            if not model.fissionableCells:
                continue
            for i, j in positions[b.parent]:
                x, y, _ = self.grid[i, j, 0].getLocalCoordinates()
                boxes.append((x - half, y - half, b.p.zbottom, x + half, y + half, b.p.ztop))
        if not boxes:
            return bounds
        lo = np.maximum(np.min(np.array(boxes)[:, :3], axis=0), bounds[:3])
        hi = np.minimum(np.max(np.array(boxes)[:, 3:], axis=0), bounds[3:])
        return tuple(lo) + tuple(hi)

    def _tallies(self, result):
        cellFilter = openmc.CellFilter(list(result.blockCells.values()))
        totals = openmc.Tally(name="block totals")
        totals.filters = [cellFilter]
        totals.scores = ["flux", "kappa-fission", "fission", "nu-fission"]
        spectrum = openmc.Tally(name="block spectrum")
        spectrum.filters = [cellFilter, openmc.EnergyFilter(result.energyBounds)]
        spectrum.scores = ["flux"]
        tallies = [totals, spectrum]
        if self.options.pinPowers:
            tallies += pinPowerTallies(result.blockModels.values())
        return openmc.Tallies(tallies)

    def _plots(self, bounds):
        x0, y0, z0, x1, y1, z1 = bounds
        width = max(x1 - x0, y1 - y0)
        zMid = 0.5 * (z0 + z1)
        xy = openmc.Plot(name="core-xy")
        xy.basis = "xy"
        xy.origin = (0.5 * (x0 + x1), 0.5 * (y0 + y1), zMid)
        xy.width = (width, width)
        xy.pixels = (1600, 1600)
        xy.color_by = "material"
        xz = openmc.Plot(name="core-xz")
        xz.basis = "xz"
        xz.origin = (0.5 * (x0 + x1), 0.5 * (y0 + y1), zMid)
        xz.width = (x1 - x0, max(z1 - z0, 1.0))
        xz.pixels = (1600, int(1600 * max(z1 - z0, 1.0) / max(x1 - x0, 1e-6)) or 400)
        xz.color_by = "material"
        return openmc.Plots([xy, xz])


def makeSettings(
    options: OpenMCOptions, bounds, runMode="eigenvalue", source=None, sourceBounds=None
):
    """OpenMC run settings: statistics, temperatures, an initial fission source and entropy."""
    s = openmc.Settings()
    s.run_mode = runMode
    s.particles = int(options.particles)
    s.batches = int(options.batches)
    if runMode == "eigenvalue":
        s.inactive = int(options.inactive)
    s.seed = int(options.seed)
    s.temperature = {
        "method": options.temperatureMethod,
        "tolerance": float(options.temperatureTolerance),
        "default": 293.6,
    }
    x0, y0, z0, x1, y1, z1 = bounds
    if source is None:
        constraints = {"fissionable": True} if runMode == "eigenvalue" else {}
        sx0, sy0, sz0, sx1, sy1, sz1 = sourceBounds or bounds
        source = openmc.IndependentSource(
            space=openmc.stats.Box((sx0, sy0, sz0), (sx1, sy1, sz1)), constraints=constraints
        )
    s.source = source
    # Thin fuel in a large box (pins, slices) rejects many sites; accept that.
    s.source_rejection_fraction = 1e-4
    if runMode == "eigenvalue":
        mesh = openmc.RegularMesh()
        mesh.lower_left = (x0, y0, z0)
        mesh.upper_right = (x1, y1, z1)
        mesh.dimension = (8, 8, 8 if z1 - z0 > 2.0 else 1)
        s.entropy_mesh = mesh
    s.output = {"tallies": False}
    # Instance offsets for every material cell are only needed for cell-instance results we never
    # ask for (distribcell tallies add their own cells). With explicit particles they are fatal:
    # each particle layer cell repeats in every compact of every block.
    s.material_cell_offsets = False
    return s


def pinPowerTallies(models):
    """A distribcell fission-power tally over the fissionable cells of every pin lattice."""
    cells = [c for m in models if m.pinLattice is not None for c in m.fissionableCells]
    tallies = []
    for cell in cells:
        t = openmc.Tally(name=f"pins {cell.id}")
        t.filters = [openmc.DistribcellFilter(cell)]
        t.scores = ["kappa-fission"]
        tallies.append(t)
    return tallies
