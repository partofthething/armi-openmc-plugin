"""
Case settings for the OpenMC plugin.

Every setting starts with ``openmc`` so the plugin's own settings are easy to find in a settings
file and in the input editor. Two framework settings pick the plugin as a solver:

* ``neutronicsKernel: OpenMC`` runs OpenMC as the global flux solver (whole-core eigenvalue, block
  power and multigroup flux mapped back onto the reactor).
* ``xsKernel: OpenMC`` runs OpenMC as the lattice physics code: one model per cross section
  group, tallied into microscopic multigroup cross sections and written as ISOTXS.
"""

from armi.physics.neutronics.settings import CONF_NEUTRONICS_KERNEL, CONF_XS_KERNEL
from armi.settings import setting

KERNEL_NAME = "OpenMC"

CONF_RUNNER = "openmcRunner"
CONF_EXECUTABLE = "openmcExecutable"
CONF_DOCKER_IMAGE = "openmcDockerImage"
CONF_CROSS_SECTIONS = "openmcCrossSections"
CONF_THREADS = "openmcThreads"
CONF_PARTICLES = "openmcParticles"
CONF_BATCHES = "openmcBatches"
CONF_INACTIVE = "openmcInactiveBatches"
CONF_SEED = "openmcSeed"
CONF_BLOCK_MODEL = "openmcBlockModel"
CONF_HETEROGENEOUS_FLAGS = "openmcHeterogeneousBlockFlags"
CONF_SYMMETRY = "openmcSymmetry"
CONF_RADIAL_BOUNDARY = "openmcRadialBoundary"
CONF_TOP_BOUNDARY = "openmcTopBoundary"
CONF_BOTTOM_BOUNDARY = "openmcBottomBoundary"
CONF_AXIAL_SLICE = "openmcAxialSlice"
CONF_TEMPERATURE_METHOD = "openmcTemperatureMethod"
CONF_TEMPERATURE_TOLERANCE = "openmcTemperatureTolerance"
CONF_THERMAL_SCATTERING = "openmcThermalScattering"
CONF_GROUP_STRUCTURE = "openmcGroupStructure"
CONF_PIN_POWERS = "openmcPinPowers"
CONF_XS_COLORSET = "openmcXSColorset"
CONF_XS_TRACE_DENSITY = "openmcXSTraceDensity"
CONF_RUN_DIR = "openmcRunDir"
CONF_PLOT = "openmcPlotGeometry"

BOUNDARIES = ["vacuum", "reflective"]


def defineSettings():
    return [
        setting.Option(KERNEL_NAME, CONF_NEUTRONICS_KERNEL),
        setting.Option(KERNEL_NAME, CONF_XS_KERNEL),
        setting.Setting(
            CONF_RUNNER,
            default="auto",
            label="How to run OpenMC",
            description=(
                "native runs the openmc executable on this machine, docker runs it in a container "
                "(the run directory and the nuclear data are mounted at the same paths), auto uses "
                "the executable when it is on the PATH and docker otherwise."
            ),
            options=["auto", "native", "docker"],
            enforcedOptions=True,
        ),
        setting.Setting(
            CONF_EXECUTABLE,
            default="openmc",
            label="OpenMC executable",
            description="Name or path of the openmc executable used by the native runner.",
        ),
        setting.Setting(
            CONF_DOCKER_IMAGE,
            default="openmc/openmc:latest",
            label="OpenMC container image",
            description="Container image used by the docker runner.",
        ),
        setting.Setting(
            CONF_CROSS_SECTIONS,
            default="",
            label="Nuclear data (cross_sections.xml)",
            description=(
                "Path to the cross_sections.xml of an OpenMC HDF5 data library. Blank uses the "
                "OPENMC_CROSS_SECTIONS environment variable."
            ),
        ),
        setting.Setting(
            CONF_THREADS,
            default=0,
            label="Threads",
            description="OpenMP threads for each OpenMC run. Zero uses every physical core.",
        ),
        setting.Setting(
            CONF_PARTICLES,
            default=10000,
            label="Particles per batch",
            description="Neutron histories per batch.",
        ),
        setting.Setting(
            CONF_BATCHES,
            default=150,
            label="Batches",
            description="Total batches, inactive ones included.",
        ),
        setting.Setting(
            CONF_INACTIVE,
            default=50,
            label="Inactive batches",
            description="Batches discarded while the fission source converges.",
        ),
        setting.Setting(
            CONF_SEED,
            default=1,
            label="Random number seed",
            description="Seed of the pseudo-random number sequence.",
        ),
        setting.Setting(
            CONF_BLOCK_MODEL,
            default="heterogeneous",
            label="Block representation",
            description=(
                "heterogeneous models every pin, duct and gap a block describes (falling back to a "
                "homogenized block when its components cannot be laid out); homogenized smears "
                "each block into one material."
            ),
            options=["heterogeneous", "homogenized"],
            enforcedOptions=True,
        ),
        setting.Setting(
            CONF_HETEROGENEOUS_FLAGS,
            default=[],
            label="Blocks to model in detail",
            description=(
                "Flags of the blocks to model heterogeneously, e.g. [fuel, control]. Empty means "
                "every block. Only used with the heterogeneous block representation."
            ),
        ),
        setting.Setting(
            CONF_SYMMETRY,
            default="native",
            label="Symmetry treatment",
            description=(
                "native keeps a quarter Cartesian core with reflective planes and unfolds other "
                "partial cores to full core; full always unfolds to full core."
            ),
            options=["native", "full"],
            enforcedOptions=True,
        ),
        setting.Setting(
            CONF_RADIAL_BOUNDARY,
            default="vacuum",
            label="Radial boundary",
            description="Boundary condition on the outer radial surface of the core model.",
            options=BOUNDARIES,
            enforcedOptions=True,
        ),
        setting.Setting(
            CONF_TOP_BOUNDARY,
            default="vacuum",
            label="Top boundary",
            description="Boundary condition at the top of the core model.",
            options=BOUNDARIES,
            enforcedOptions=True,
        ),
        setting.Setting(
            CONF_BOTTOM_BOUNDARY,
            default="vacuum",
            label="Bottom boundary",
            description="Boundary condition at the bottom of the core model.",
            options=BOUNDARIES,
            enforcedOptions=True,
        ),
        setting.Setting(
            CONF_AXIAL_SLICE,
            default=-1.0,
            label="2D slice elevation (cm)",
            description=(
                "Negative builds the full 3D core. A non-negative elevation builds a 2D model of "
                "the blocks at that height, with reflective top and bottom."
            ),
        ),
        setting.Setting(
            CONF_TEMPERATURE_METHOD,
            default="interpolation",
            label="Temperature treatment",
            description="How OpenMC uses data at temperatures between the library's.",
            options=["interpolation", "nearest"],
            enforcedOptions=True,
        ),
        setting.Setting(
            CONF_TEMPERATURE_TOLERANCE,
            default=50.0,
            label="Temperature tolerance (K)",
            description=(
                "How far a material temperature may sit outside the library's temperatures, or "
                "from the nearest one, and still use it."
            ),
        ),
        setting.Setting(
            CONF_THERMAL_SCATTERING,
            default="auto",
            label="Thermal scattering (S(a,b))",
            description=(
                "auto adds S(a,b) tables to materials recognized as water, heavy water, graphite, "
                "beryllium, beryllium oxide or zirconium hydride; none never adds them."
            ),
            options=["auto", "none"],
            enforcedOptions=True,
        ),
        setting.Setting(
            CONF_GROUP_STRUCTURE,
            default="",
            label="Tally group structure",
            description=(
                "Energy groups of the multigroup flux and cross section tallies. Blank uses the "
                "case's groupStructure setting."
            ),
        ),
        setting.Setting(
            CONF_PIN_POWERS,
            default=False,
            label="Tally pin powers",
            description="Tally the fission power of every pin and store it on the blocks.",
        ),
        setting.Setting(
            CONF_XS_COLORSET,
            default=True,
            label="Use driver colorsets",
            description=(
                "When a cross section group names a driverID, surround the representative block "
                "with driver blocks instead of modeling an infinite lattice of it alone."
            ),
        ),
        setting.Setting(
            CONF_XS_TRACE_DENSITY,
            default=1e-12,
            label="Trace density (atoms/b-cm)",
            description=(
                "Density at which nuclides missing from a block are added so that cross sections "
                "can still be tallied for them."
            ),
        ),
        setting.Setting(
            CONF_RUN_DIR,
            default="openmc-runs",
            label="Run directory",
            description="Directory, relative to the case, where each OpenMC run gets a folder.",
        ),
        setting.Setting(
            CONF_PLOT,
            default=False,
            label="Plot geometry",
            description="Write PNG slices of each model geometry beside its inputs.",
        ),
    ]


def defineSettingEditors(cs):
    """Group the plugin's settings into one section of the input editor."""
    group = "OpenMC"
    hint = (
        "OpenMC Monte Carlo, used as the global flux solver (neutronicsKernel: OpenMC) or as the "
        "lattice physics code that makes ISOTXS cross sections (xsKernel: OpenMC)."
    )
    common = {"group": group}
    editors = {
        CONF_RUNNER: {**common, "order": 1, "groupHint": hint},
        CONF_EXECUTABLE: {
            **common,
            "order": 2,
            "widget": "path",
            "showWhen": {"setting": CONF_RUNNER, "in": ["auto", "native"]},
        },
        CONF_DOCKER_IMAGE: {
            **common,
            "order": 2,
            "showWhen": {"setting": CONF_RUNNER, "in": ["auto", "docker"]},
        },
        CONF_CROSS_SECTIONS: {**common, "order": 3, "widget": "path"},
        CONF_THREADS: {**common, "order": 4},
        CONF_PARTICLES: {**common, "order": 10},
        CONF_BATCHES: {**common, "order": 11},
        CONF_INACTIVE: {**common, "order": 12},
        CONF_SEED: {**common, "order": 13},
        CONF_BLOCK_MODEL: {**common, "order": 20},
        CONF_HETEROGENEOUS_FLAGS: {
            **common,
            "order": 21,
            "widget": "choices",
            "options": ["fuel", "control", "shield", "reflector", "plenum", "duct"],
            "allowOther": True,
            "showWhen": {"setting": CONF_BLOCK_MODEL, "equals": "heterogeneous"},
        },
        CONF_SYMMETRY: {**common, "order": 30},
        CONF_RADIAL_BOUNDARY: {**common, "order": 31},
        CONF_TOP_BOUNDARY: {**common, "order": 32},
        CONF_BOTTOM_BOUNDARY: {**common, "order": 33},
        CONF_AXIAL_SLICE: {**common, "order": 34},
        CONF_TEMPERATURE_METHOD: {**common, "order": 40},
        CONF_TEMPERATURE_TOLERANCE: {**common, "order": 41},
        CONF_THERMAL_SCATTERING: {**common, "order": 42},
        CONF_GROUP_STRUCTURE: {**common, "order": 50},
        CONF_PIN_POWERS: {**common, "order": 51},
        CONF_XS_COLORSET: {**common, "order": 60},
        CONF_XS_TRACE_DENSITY: {**common, "order": 61},
        CONF_RUN_DIR: {**common, "order": 70, "widget": "path", "pathKind": "dir"},
        CONF_PLOT: {**common, "order": 71},
    }
    return editors
