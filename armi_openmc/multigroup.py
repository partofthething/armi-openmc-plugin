"""
OpenMC in multigroup mode, on the cross sections of an ARMI ISOTXS library.

The continuous-energy model answers "what does the reactor do". The multigroup model answers a
narrower question: what does the *multigroup problem* a deterministic core solver is given do, when
solved without its approximations (diffusion, SP\\ :sub:`N`, coarse mesh)? Every block is smeared
and gets the macroscopic cross sections of its own number densities times the microscopic cross
sections of its cross section ID (``<nuclide label><xsID>`` in the library), exactly what an
ISOTXS-fed nodal code builds. A difference between this run and the deterministic core is the core
solver's; a difference between this run and the continuous-energy one is the cross sections'.

Per block composition, :py:class:`MultigroupMaterialFactory` makes one ``openmc.XSdata``:

* total and absorption (fission plus every capture reaction), flux-weighted as in the ISOTXS;
* the P\\ :sub:`0` and P\\ :sub:`1` scattering matrices (Legendre order 1), as *production*
  matrices -- the ISOTXS total scattering matrix counts the two neutrons of (n,2n) -- with a
  multiplicity matrix that makes total = absorption + scattering;
* ν-fission as a full group-to-group matrix, the sum over nuclides of N ν σ\\ :sub:`f` χ, so
  that each nuclide keeps its own fission spectrum;
* fission and κ-fission (``kappa`` per fission from the library, or 200 MeV where it gives none,
  as DRAGON's ISOTXS does), for the power tallies.
"""

import hashlib

import numpy as np
import openmc
import openmc.mgxs
from armi import runLog
from armi.nucDirectory import nuclideBases

#: Energy per fission (eV) where the library carries none.
DEFAULT_KAPPA_EV = 200.0e6
#: Lower bound (eV) of the lowest group: an ISOTXS gives only upper bounds.
LOWEST_ENERGY_EV = 1.0e-5
#: Temperature (K) of the multigroup data: the default of ``makeSettings``, so every material
#: (which has no temperature of its own) picks it.
DATA_TEMPERATURE = 293.6
CAPTURES = ("nGamma", "nalph", "np", "nd", "nt")


def energyBounds(lib) -> np.ndarray:
    """
    Group bounds of an ARMI library, ascending, in eV.

    A library condensed onto a structure whose bounds do not all fall between its fine groups
    (ANL116 from SHEM-361, say) has empty groups of zero width. OpenMC needs strictly increasing
    bounds, so each repeated bound is moved up by a part in a billion. Nothing reaches those
    groups anyway: their cross sections and their χ are zero.
    """
    upper = np.asarray(lib.neutronEnergyUpperBounds, dtype=float)
    bounds = np.concatenate([[LOWEST_ENERGY_EV], upper[::-1]])
    repeated = 0
    for i in range(1, len(bounds)):
        if bounds[i] <= bounds[i - 1]:
            bounds[i] = bounds[i - 1] * (1.0 + 1.0e-9)
            repeated += 1
    if repeated:
        runLog.warning(
            f"OpenMC multigroup: the library has {repeated} empty (zero-width) groups; "
            "their bounds are moved apart by a part in 1e9",
            single=True,
        )
    return bounds


def _vector(values, nGroups):
    if values is None:
        return np.zeros(nGroups)
    # an empty group of a library can hold NaN (zero flux weight): it carries nothing
    return np.nan_to_num(np.asarray(values, dtype=float).reshape(nGroups, -1)[:, 0])


def _matrix(values, nGroups):
    """An ARMI scattering matrix ([to, from]) as a dense array, or zeros."""
    if values is None:
        return np.zeros((nGroups, nGroups))
    dense = values.toarray() if hasattr(values, "toarray") else np.asarray(values, dtype=float)
    return np.nan_to_num(dense.reshape(nGroups, nGroups))


class MultigroupMaterialFactory:
    """
    OpenMC materials of whole blocks on macroscopic multigroup data from an ARMI library.

    Parameters
    ----------
    lib : armi.nuclearDataIO.xsLibraries.IsotxsLibrary
        The microscopic cross sections, named ``<nuclide label><xsID>``.
    traceDensity : float
        Number densities below this (atoms/barn-cm) are left out.
    """

    def __init__(self, lib, traceDensity=0.0):
        self.lib = lib
        self.nGroups = lib.numGroups
        self.bounds = energyBounds(lib)
        self.groups = openmc.mgxs.EnergyGroups(self.bounds)
        self.library = openmc.MGXSLibrary(self.groups)
        self.traceDensity = traceDensity
        self.materials = []
        self._xsdata = {}
        self._fissionable = set()
        self._missing = {}
        self._digest = hashlib.sha1()

    # ----------------------------------------------------------------------------------------
    def forBlock(self, b) -> openmc.Material:
        """A material for block ``b``: its own, on the data of its composition (shared)."""
        suffix = b.getMicroSuffix()
        densities = {n: d for n, d in b.getNumberDensities().items() if d > self.traceDensity}
        key = (suffix, tuple(sorted((n, round(d, 12)) for n, d in densities.items())))
        if key not in self._xsdata:
            self._xsdata[key] = self._macroscopic(
                f"{suffix}-{len(self._xsdata) + 1}", suffix, densities
            )
        name = self._xsdata[key]
        mat = openmc.Material(name=b.getName())
        mat.set_density("macro", 1.0)
        mat.add_macroscopic(name)
        if name in self._fissionable:
            self._fissionable.add(mat.id)
        self.materials.append(mat)
        return mat

    def isFissionable(self, mat) -> bool:
        return mat is not None and mat.id in self._fissionable

    def make(self, name, numberDensities, temperatureInC, reuse=False):
        raise ValueError(
            f"OpenMC multigroup mode smears whole blocks; {name} would need a component material "
            "(use openmcBlockModel: homogenized)"
        )

    @property
    def digest(self) -> str:
        """A fingerprint of every cross section written, for reusing a finished statepoint."""
        return self._digest.hexdigest()

    def report(self):
        runLog.info(
            f"OpenMC multigroup: {len(self.materials)} block materials on {len(self._xsdata)} "
            f"macroscopic cross section sets, {self.nGroups} groups"
        )
        for label, density in sorted(self._missing.items(), key=lambda kv: -kv[1]):
            runLog.warning(
                f"OpenMC multigroup: {label} is not in the library and is left out "
                f"(up to {density:.3e} atoms/barn-cm)",
                single=True,
                label=f"openmc-mg-missing-{label}",
            )

    def export(self, path):
        self.library.export_to_hdf5(str(path))

    # ----------------------------------------------------------------------------------------
    def _macroscopic(self, name, suffix, densities) -> str:
        G = self.nGroups
        total, absorption, fission, kappa = (np.zeros(G) for _ in range(4))
        nuFission = np.zeros((G, G))  # [from, to]
        scatter = np.zeros((2, G, G))  # [order, to, from]
        for nucName, density in densities.items():
            base = nuclideBases.byName.get(nucName)
            label = f"{base.label if base is not None else nucName}{suffix}"
            if label not in self.lib:
                self._missing[label] = max(density, self._missing.get(label, 0.0))
                continue
            nuclide = self.lib[label]
            micros = nuclide.micros
            total += density * _vector(micros.total, G)
            capture = sum(_vector(getattr(micros, r, None), G) for r in CAPTURES)
            sigF = _vector(micros.fission, G)
            absorption += density * (capture + sigF)
            if sigF.any():
                nu = _vector(micros.neutronsPerFission, G)
                chi = np.asarray(micros.chi, dtype=float) if micros.chi is not None else np.zeros(G)
                if chi.ndim == 2:  # a fission matrix [to, from]
                    nuFission += density * (chi * (nu * sigF)[None, :]).T
                else:
                    nuFission += density * np.outer(nu * sigF, chi.reshape(G))
                fission += density * sigF
                try:
                    ev = float(nuclide.isotxsMetadata["efiss"] or 0.0)
                except (KeyError, TypeError):
                    ev = 0.0
                # an ISOTXS gives watt-seconds per fission
                kappaEV = ev / 1.602176634e-19 if ev > 0.0 else DEFAULT_KAPPA_EV
                kappa += density * sigF * kappaEV
            scatter[0] += density * _matrix(micros.totalScatter, G)
            higher = micros.higherOrderScatter or {}
            if 1 in higher:
                scatter[1] += density * _matrix(higher[1], G)

        # The production matrix counts (n,2n) twice; the multiplicity matrix makes the collision
        # probabilities right (total = absorption + scattering).
        production = scatter[0].sum(axis=0)  # out of each group
        scattering = np.clip(total - absorption, 0.0, None)
        multiplicity = np.ones(G)
        nonzero = scattering > 0.0
        multiplicity[nonzero] = production[nonzero] / scattering[nonzero]
        multiplicity = np.where(multiplicity > 0.0, multiplicity, 1.0)

        xs = openmc.XSdata(name, self.groups, temperatures=[DATA_TEMPERATURE])
        xs.order = 1
        xs.set_total(total, temperature=DATA_TEMPERATURE)
        xs.set_absorption(absorption, temperature=DATA_TEMPERATURE)
        # OpenMC wants [from, to, order]
        xs.set_scatter_matrix(np.transpose(scatter, (2, 1, 0)), temperature=DATA_TEMPERATURE)
        xs.set_multiplicity_matrix(
            np.repeat(multiplicity[:, None], G, axis=1), temperature=DATA_TEMPERATURE
        )
        if fission.any():
            xs.set_fission(fission, temperature=DATA_TEMPERATURE)
            xs.set_nu_fission(nuFission, temperature=DATA_TEMPERATURE)
            xs.set_kappa_fission(kappa, temperature=DATA_TEMPERATURE)
            self._fissionable.add(name)
        self.library.add_xsdata(xs)
        for array in (total, absorption, scatter, nuFission, kappa):
            self._digest.update(np.ascontiguousarray(array).tobytes())
        return name
