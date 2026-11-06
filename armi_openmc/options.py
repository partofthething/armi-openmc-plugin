"""Modeling and execution options, gathered from case settings or set directly in code."""

import os
from dataclasses import dataclass, field, fields

from armi_openmc import settings as S


@dataclass
class OpenMCOptions:
    """Everything that shapes an OpenMC model and its run, independent of the reactor."""

    particles: int = 10000
    batches: int = 150
    inactive: int = 50
    seed: int = 1
    threads: int = 0

    runner: str = "auto"
    executable: str = "openmc"
    dockerImage: str = "openmc/openmc:latest"
    crossSections: str = ""

    blockModel: str = "heterogeneous"
    heterogeneousFlags: list[str] = field(default_factory=list)
    symmetry: str = "native"
    radialBoundary: str = "vacuum"
    topBoundary: str = "vacuum"
    bottomBoundary: str = "vacuum"
    axialSlice: float | None = None
    coreRadius: float = 0.0
    outerFill: str = ""

    temperatureMethod: str = "interpolation"
    temperatureTolerance: float = 50.0
    thermalScattering: str = "auto"

    groupStructure: str = "ANL33"
    pinPowers: bool = False

    particleModel: str = "explicit"
    particleSegmentLength: float = 5.0
    particleSeed: int = 1
    particleMeshPitch: float = 0.0

    xsColorset: bool = True
    xsTraceDensity: float = 1e-12

    runDir: str = "openmc-runs"
    plotGeometry: bool = False
    reuseStatepoint: bool = False

    energyMode: str = "continuous-energy"
    multigroupLibrary: str = ""

    _SETTINGS = {
        "particles": S.CONF_PARTICLES,
        "batches": S.CONF_BATCHES,
        "inactive": S.CONF_INACTIVE,
        "seed": S.CONF_SEED,
        "threads": S.CONF_THREADS,
        "runner": S.CONF_RUNNER,
        "executable": S.CONF_EXECUTABLE,
        "dockerImage": S.CONF_DOCKER_IMAGE,
        "crossSections": S.CONF_CROSS_SECTIONS,
        "blockModel": S.CONF_BLOCK_MODEL,
        "heterogeneousFlags": S.CONF_HETEROGENEOUS_FLAGS,
        "symmetry": S.CONF_SYMMETRY,
        "radialBoundary": S.CONF_RADIAL_BOUNDARY,
        "topBoundary": S.CONF_TOP_BOUNDARY,
        "bottomBoundary": S.CONF_BOTTOM_BOUNDARY,
        "axialSlice": S.CONF_AXIAL_SLICE,
        "coreRadius": S.CONF_CORE_RADIUS,
        "outerFill": S.CONF_OUTER_FILL,
        "temperatureMethod": S.CONF_TEMPERATURE_METHOD,
        "temperatureTolerance": S.CONF_TEMPERATURE_TOLERANCE,
        "thermalScattering": S.CONF_THERMAL_SCATTERING,
        "groupStructure": S.CONF_GROUP_STRUCTURE,
        "pinPowers": S.CONF_PIN_POWERS,
        "particleModel": S.CONF_PARTICLE_MODEL,
        "particleSegmentLength": S.CONF_PARTICLE_SEGMENT,
        "particleSeed": S.CONF_PARTICLE_SEED,
        "particleMeshPitch": S.CONF_PARTICLE_MESH_PITCH,
        "xsColorset": S.CONF_XS_COLORSET,
        "xsTraceDensity": S.CONF_XS_TRACE_DENSITY,
        "runDir": S.CONF_RUN_DIR,
        "plotGeometry": S.CONF_PLOT,
        "reuseStatepoint": S.CONF_REUSE,
        "energyMode": S.CONF_ENERGY_MODE,
        "multigroupLibrary": S.CONF_MG_LIBRARY,
    }

    @classmethod
    def fromSettings(cls, cs, **overrides):
        """Build options from case settings; keyword arguments win over the settings."""
        values = {}
        for attr, name in cls._SETTINGS.items():
            if name in cs:
                values[attr] = cs[name]
        if not values.get("groupStructure") and "groupStructure" in cs:
            values["groupStructure"] = cs["groupStructure"]
        slice_ = values.get("axialSlice")
        values["axialSlice"] = None if slice_ is None or slice_ < 0 else float(slice_)
        values["heterogeneousFlags"] = list(values.get("heterogeneousFlags") or [])
        if values.get("multigroupLibrary") and not os.path.isabs(values["multigroupLibrary"]):
            values["multigroupLibrary"] = os.path.join(
                cs.inputDirectory, values["multigroupLibrary"]
            )
        if values.get("runDir") and not os.path.isabs(values["runDir"]):
            values["runDir"] = os.path.join(cs.inputDirectory, values["runDir"])
        values.update(overrides)
        known = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in values.items() if k in known})

    @property
    def multigroup(self) -> bool:
        return self.energyMode == "multigroup"

    @property
    def crossSectionsPath(self) -> str:
        """The cross_sections.xml this run uses, from the option, the environment or the default."""
        path = self.crossSections or os.environ.get("OPENMC_CROSS_SECTIONS", "")
        if not path or not os.path.exists(path):
            raise FileNotFoundError(
                f"OpenMC nuclear data not found at {path}. Set the {S.CONF_CROSS_SECTIONS} setting "
                "or the OPENMC_CROSS_SECTIONS environment variable."
            )
        return os.path.abspath(path)
