"""
Turn ARMI number densities into OpenMC materials.

ARMI names nuclides its own way (``U235``, ``AM242M``, elemental ``FE``, lumped ``LFP35``) and
holds number densities in atoms/barn-cm. OpenMC wants GNDS names (``U235``, ``Am242_m1``) that
exist in the nuclear data library. :py:class:`NuclideMapper` does the translation, expanding
natural elements over the isotopes the library has and dropping, with a warning, what the
library cannot represent.
"""

import functools
import math
from collections import defaultdict

import openmc
import openmc.data
from armi import runLog
from armi.nucDirectory import nuclideBases

KELVIN_OFFSET = 273.15


@functools.lru_cache(maxsize=4)
def _libraryContents(crossSectionsPath):
    """Neutron nuclide names and thermal scattering table names in a cross_sections.xml."""
    library = openmc.data.DataLibrary.from_xml(crossSectionsPath)
    neutron, thermal, paths = set(), set(), {}
    for entry in library.libraries:
        if entry["type"] == "neutron":
            neutron.update(entry["materials"])
            for name in entry["materials"]:
                paths[name] = entry["path"]
        elif entry["type"] == "thermal":
            thermal.update(entry["materials"])
    return frozenset(neutron), frozenset(thermal), paths


class NuclideMapper:
    """Maps ARMI nuclide names onto the nuclides of one OpenMC data library."""

    def __init__(self, crossSectionsPath):
        self.crossSectionsPath = crossSectionsPath
        self.available, self.thermalTables, self.dataPaths = _libraryContents(crossSectionsPath)
        self._cache = {}
        self._warned = set()

    def expand(self, armiName):
        """
        OpenMC nuclides and fractions that make up one ARMI nuclide.

        Returns an empty list for nuclides OpenMC cannot represent (lumped fission products,
        dummy nuclides, isotopes missing from the library).
        """
        if armiName not in self._cache:
            self._cache[armiName] = self._expand(armiName)
        return self._cache[armiName]

    def _expand(self, armiName):
        base = nuclideBases.byName.get(armiName)
        if base is None or isinstance(
            base, (nuclideBases.LumpNuclideBase, nuclideBases.DummyNuclideBase)
        ):
            self._warn(armiName, "is a lumped or dummy nuclide with no continuous-energy data")
            return []
        symbol = openmc.data.ATOMIC_SYMBOL[base.z]
        if isinstance(base, nuclideBases.NaturalNuclideBase):
            return self._expandElement(armiName, symbol)
        name = f"{symbol}{base.a}" + (f"_m{base.state}" if base.state else "")
        if name in self.available:
            return [(name, 1.0)]
        self._warn(armiName, f"({name}) is not in the nuclear data library")
        return []

    def _expandElement(self, armiName, symbol):
        abundances = {
            nuc: frac
            for nuc, frac in openmc.data.NATURAL_ABUNDANCE.items()
            if openmc.data.zam(nuc)[0] == openmc.data.ATOMIC_NUMBER[symbol]
        }
        present = {nuc: frac for nuc, frac in abundances.items() if nuc in self.available}
        if not present:
            self._warn(armiName, "has no natural isotopes in the nuclear data library")
            return []
        missing = sum(abundances.values()) - sum(present.values())
        if missing > 1e-6:
            self._warn(
                armiName,
                f"is missing isotopes worth {missing:.2%} of the element; "
                "the rest are renormalized",
            )
        total = sum(present.values())
        return [(nuc, frac / total) for nuc, frac in sorted(present.items())]

    def _warn(self, armiName, why):
        if armiName not in self._warned:
            self._warned.add(armiName)
            runLog.warning(f"OpenMC: nuclide {armiName} {why}, so it is left out of the model.")

    def toOpenMC(self, numberDensities):
        """ARMI {name: atoms/b-cm} to OpenMC {name: atoms/b-cm}, summing shared isotopes."""
        out = defaultdict(float)
        for armiName, density in numberDensities.items():
            if density <= 0.0:
                continue
            for name, frac in self.expand(armiName):
                out[name] += density * frac
        return dict(out)


# Thermal scattering tables, keyed by a test on the atom fractions of a material's elements.
def _thermalTables(elementFractions, available):
    """Pick S(a,b) tables for a material from what it is made of."""
    f = elementFractions
    h, o, c, be, zr, d = (f.get(e, 0.0) for e in ("H", "O", "C", "Be", "Zr", "H2"))
    tables = []
    if h > 0.0 and o > 0.0 and h / o > 1.0 and zr < h / 4:
        tables.append("c_H_in_H2O")
    elif h > 0.0 and zr > 0.0 and h / zr > 1.0:
        tables += ["c_H_in_ZrH", "c_Zr_in_ZrH"]
    if d > 0.0 and o > 0.0 and d / o > 1.0:
        tables.append("c_D_in_D2O")
    if c > 0.9:
        tables.append("c_Graphite")
    if be > 0.9:
        tables.append("c_Be")
    elif be > 0.3 and o > 0.3 and abs(be - o) < 0.1:
        tables += ["c_Be_in_BeO", "c_O_in_BeO"]
    return [t for t in tables if t in available]


class MaterialFactory:
    """
    Builds OpenMC materials from ARMI number densities, one call per region.

    Materials are not shared between ARMI objects: every component of every block gets its own,
    so that later depletion or temperature coupling can change one without touching the rest.
    """

    def __init__(self, mapper: NuclideMapper, thermalScattering="auto", traceNuclides=None):
        self.mapper = mapper
        self.thermalScattering = thermalScattering
        #: {ARMI name: density} added to every material that lacks the nuclide (XS generation).
        self.traceNuclides = traceNuclides or {}
        self.materials = []
        self._shared = {}

    def make(self, name, numberDensities, temperatureInC, reuse=False):
        """
        A material for these ARMI number densities, or None for a void.

        Parameters
        ----------
        name : str
            Label for the material (shows up in OpenMC output and plots).
        numberDensities : dict
            ARMI nuclide name to atoms/barn-cm.
        temperatureInC : float
            Material temperature in Celsius.
        reuse : bool
            Return the material already made for identical densities and temperature, if any.
            For regions that are deliberately shared (particle layers of identical compacts),
            never for ones that may later need to change on their own.
        """
        densities = {k: v for k, v in numberDensities.items() if v > 0.0}
        for nuc, trace in self.traceNuclides.items():
            densities.setdefault(nuc, trace)
        omc = self.mapper.toOpenMC(densities)
        if not omc or sum(omc.values()) <= 0.0:
            return None
        key = None
        if reuse:
            key = (
                tuple(sorted((n, round(d, 12)) for n, d in omc.items())),
                round(temperatureInC, 3),
            )
            if key in self._shared:
                return self._shared[key]
        mat = openmc.Material(name=name, temperature=temperatureInC + KELVIN_OFFSET)
        for nuc, density in sorted(omc.items()):
            mat.add_nuclide(nuc, density, "ao")
        mat.set_density("sum")
        if self.thermalScattering == "auto":
            for table in _thermalTables(_elementFractions(omc), self.mapper.thermalTables):
                mat.add_s_alpha_beta(table)
        self.materials.append(mat)
        if key is not None:
            self._shared[key] = mat
        return mat


def _elementFractions(omcDensities):
    """Atom fraction of each element, with deuterium counted apart from hydrogen as ``H2``."""
    total = sum(omcDensities.values())
    out = defaultdict(float)
    for nuc, dens in omcDensities.items():
        z, a, _m = openmc.data.zam(nuc)
        key = "H2" if (z == 1 and a == 2) else openmc.data.ATOMIC_SYMBOL[z]
        out[key] += dens / total
    return out


def mixNumberDensities(parts):
    """
    Volume-weighted mix of several regions' number densities.

    Parameters
    ----------
    parts : list of (dict, float, float)
        (number densities, area or volume, temperature in C) for each region.

    Returns
    -------
    (dict, float)
        The mixed number densities and the volume-weighted temperature.
    """
    total = sum(v for _, v, _ in parts)
    if total <= 0.0:
        return {}, (parts[0][2] if parts else 20.0)
    mixed = defaultdict(float)
    for ndens, vol, _ in parts:
        for nuc, dens in ndens.items():
            mixed[nuc] += dens * vol / total
    temp = sum(t * v for _, v, t in parts) / total
    return dict(mixed), temp


def isClose(a, b, rel=1e-6, abs_=1e-6):
    return math.isclose(a, b, rel_tol=rel, abs_tol=abs_)
