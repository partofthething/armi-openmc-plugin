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
CONF_CORE_RADIUS = "openmcCoreRadius"
CONF_OUTER_FILL = "openmcOuterFill"
CONF_TEMPERATURE_METHOD = "openmcTemperatureMethod"
CONF_TEMPERATURE_TOLERANCE = "openmcTemperatureTolerance"
CONF_THERMAL_SCATTERING = "openmcThermalScattering"
CONF_GROUP_STRUCTURE = "openmcGroupStructure"
CONF_PIN_POWERS = "openmcPinPowers"
CONF_XS_COLORSET = "openmcXSColorset"
CONF_XS_TRACE_DENSITY = "openmcXSTraceDensity"
CONF_PARTICLE_MODEL = "openmcParticleModel"
CONF_PARTICLE_SEGMENT = "openmcParticleSegmentLength"
CONF_PARTICLE_SEED = "openmcParticleSeed"
CONF_PARTICLE_MESH_PITCH = "openmcParticleMeshPitch"
CONF_RUN_DIR = "openmcRunDir"
CONF_PLOT = "openmcPlotGeometry"
CONF_REUSE = "openmcReuseStatepoint"
CONF_ENERGY_MODE = "openmcEnergyMode"
CONF_MG_LIBRARY = "openmcMultigroupLibrary"
CONF_MG_SCATTERING = "openmcMultigroupScattering"
CONF_DEPLETION = "openmcDepletion"
CONF_DEPLETION_CHAIN = "openmcDepletionChain"
CONF_DEPLETION_INTEGRATOR = "openmcDepletionIntegrator"
CONF_PYTHON = "openmcPython"
CONF_DEPLETION_RATES = "openmcDepletionReactionRates"
CONF_DEPLETION_GROUPS = "openmcDepletionGroups"

ENERGY_MODES = ["continuous-energy", "multigroup"]
MG_SCATTERING = ["legendre", "transport-corrected"]
#: openmc.deplete integrators by setting value; cecm (CE/CM) is OpenMC's usual choice.
INTEGRATORS = {
    "predictor": "PredictorIntegrator",
    "cecm": "CECMIntegrator",
    "celi": "CELIIntegrator",
    "leqi": "LEQIIntegrator",
    "epc-rk4": "EPCRK4Integrator",
    "si-celi": "SICELIIntegrator",
}

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
            CONF_CORE_RADIUS,
            default=0.0,
            label="Core radius (cm)",
            description=(
                "Radius of the cylinder that bounds the core model. Zero: a hexagonal core is "
                "bounded just around its outermost assemblies, a Cartesian one by the box of its "
                "lattice. A Cartesian core inside a larger cylinder (a pool, a drum) has the "
                "outer fill between its assemblies and the cylinder."
            ),
        ),
        setting.Setting(
            CONF_OUTER_FILL,
            default="",
            label="Fill outside the assemblies",
            description=(
                "Block type whose homogenized composition fills the space between the assembly "
                "lattice and the core boundary, e.g. a permanent side reflector or the water of a "
                "pool. Taken from the core, or else from an assembly design in the blueprints "
                "that is not loaded. Blank: void."
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
            CONF_PARTICLE_MODEL,
            default="explicit",
            label="Particle fuel",
            description=(
                "explicit packs every TRISO-like particle (from the particles section of the "
                "blueprints) randomly in its host on a Cartesian mesh; homogenized smears matrix "
                "and particles together. Only used for heterogeneous blocks."
            ),
            options=["explicit", "homogenized"],
            enforcedOptions=True,
        ),
        setting.Setting(
            CONF_PARTICLE_SEGMENT,
            default=5.0,
            label="Particle packing segment (cm)",
            description=(
                "Hosts are packed in segments of about this length, stacked up the block: like "
                "compacts stacked in a fuel hole, and much faster to pack than a whole block."
            ),
        ),
        setting.Setting(
            CONF_PARTICLE_SEED,
            default=1,
            label="Particle packing seed",
            description="Random seed of the particle packing.",
        ),
        setting.Setting(
            CONF_PARTICLE_MESH_PITCH,
            default=0.0,
            label="Particle mesh pitch (cm)",
            description=(
                "Pitch of the Cartesian mesh the particles are sorted onto for fast tracking. "
                "Zero picks one with a few particles per mesh element."
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
        setting.Setting(
            CONF_REUSE,
            default=False,
            label="Reuse a finished statepoint",
            description=(
                "Read the statepoint already in the run directory instead of running OpenMC, "
                "when its model.xml is identical to the one this run builds: recovers a run "
                "whose OpenMC finished but whose ARMI side did not."
            ),
        ),
        setting.Setting(
            CONF_ENERGY_MODE,
            default="continuous-energy",
            label="Energy treatment",
            description=(
                "continuous-energy: OpenMC's own nuclear data (openmcCrossSections). multigroup: "
                "every block smeared, on macroscopic cross sections built from an ARMI ISOTXS "
                "library (the one the case's lattice physics made, or openmcMultigroupLibrary) "
                "and the block's number densities; solves exactly the multigroup problem a "
                "deterministic core solver is given, without its spatial or angular approximations."
            ),
            options=ENERGY_MODES,
        ),
        setting.Setting(
            CONF_MG_LIBRARY,
            default="",
            label="Multigroup library",
            description=(
                "ISOTXS file for the multigroup mode, relative to the case. Empty: the library "
                "on the core (made by the case's cross section interface)."
            ),
        ),
        setting.Setting(
            CONF_MG_SCATTERING,
            default="legendre",
            label="Multigroup scattering",
            description=(
                "legendre: every Legendre order of the ISOTXS. OpenMC clips the negative part of "
                "each group-to-group kernel, which in fine group structures (hundreds of groups "
                "and up) makes forward scattering far too isotropic and the core leak too little. "
                "transport-corrected: isotropic scattering with each group's P1 out-scattering "
                "removed from its total and its self-scatter (outflow transport correction); "
                "positive by construction, and the better choice for fine groups."
            ),
            options=MG_SCATTERING,
        ),
        setting.Setting(
            CONF_DEPLETION,
            default=False,
            label="Deplete with OpenMC",
            description=(
                "Burn the reactor with openmc.deplete: at every time node, the depletable "
                "components are depleted over the step that just ended (transport and Bateman "
                "solution in one OpenMC run) and their number densities, every nuclide of the "
                "chain included, are written back onto them. Downtime between cycles is decayed "
                "at zero power. Independent of the global flux kernel."
            ),
        ),
        setting.Setting(
            CONF_DEPLETION_CHAIN,
            default="",
            label="Depletion chain",
            description=(
                "OpenMC depletion chain XML, matched to the nuclear data library. Blank uses the "
                "OPENMC_CHAIN_FILE environment variable."
            ),
        ),
        setting.Setting(
            CONF_DEPLETION_INTEGRATOR,
            default="cecm",
            label="Depletion integrator",
            description=(
                "openmc.deplete time integrator: predictor (one transport per step), cecm "
                "(CE/CM predictor-corrector, two), celi, leqi, epc-rk4 or si-celi."
            ),
            options=list(INTEGRATORS),
            enforcedOptions=True,
        ),
        setting.Setting(
            CONF_DEPLETION_RATES,
            default="flux",
            label="Depletion reaction rates",
            description=(
                "flux: tally the flux in openmcDepletionGroups and collapse each nuclide's cross "
                "sections with it, with direct tallies kept for fission and capture of the major "
                "actinides (where self-shielding matters most); many times faster in burned fuel, "
                "which carries hundreds of nuclides. direct: tally every reaction of every nuclide."
            ),
            options=["flux", "direct"],
            enforcedOptions=True,
        ),
        setting.Setting(
            CONF_DEPLETION_GROUPS,
            default="SHEM-361",
            label="Depletion flux groups",
            description=(
                "Group structure of the flux the reaction rates are collapsed with, by its OpenMC "
                "name (openmc.mgxs.GROUP_STRUCTURES), e.g. SHEM-361, CASMO-70, VITAMIN-J-175."
            ),
        ),
        setting.Setting(
            CONF_PYTHON,
            default="python",
            label="Python with OpenMC",
            description=(
                "Python interpreter with the openmc package and its shared library, used by the "
                "native runner for depletion. The docker runner uses the image's own."
            ),
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
    particleOnly = {"showWhen": {"setting": CONF_PARTICLE_MODEL, "equals": "explicit"}}
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
        CONF_CORE_RADIUS: {**common, "order": 35},
        CONF_OUTER_FILL: {**common, "order": 36},
        CONF_TEMPERATURE_METHOD: {**common, "order": 40},
        CONF_TEMPERATURE_TOLERANCE: {**common, "order": 41},
        CONF_THERMAL_SCATTERING: {**common, "order": 42},
        CONF_GROUP_STRUCTURE: {**common, "order": 50},
        CONF_PIN_POWERS: {**common, "order": 51},
        CONF_XS_COLORSET: {**common, "order": 60},
        CONF_XS_TRACE_DENSITY: {**common, "order": 61},
        CONF_PARTICLE_MODEL: {**common, "order": 25},
        CONF_PARTICLE_SEGMENT: {**common, "order": 26, **particleOnly},
        CONF_PARTICLE_SEED: {**common, "order": 27, **particleOnly},
        CONF_PARTICLE_MESH_PITCH: {**common, "order": 28, **particleOnly},
        CONF_RUN_DIR: {**common, "order": 70, "widget": "path", "pathKind": "dir"},
        CONF_PLOT: {**common, "order": 71},
        CONF_ENERGY_MODE: {**common, "order": 5},
        CONF_MG_LIBRARY: {
            **common,
            "order": 6,
            "widget": "path",
            "showWhen": {"setting": CONF_ENERGY_MODE, "equals": "multigroup"},
        },
        CONF_MG_SCATTERING: {
            **common,
            "order": 7,
            "showWhen": {"setting": CONF_ENERGY_MODE, "equals": "multigroup"},
        },
        CONF_DEPLETION: {**common, "order": 80},
        CONF_DEPLETION_CHAIN: {
            **common,
            "order": 81,
            "widget": "path",
            "showWhen": {"setting": CONF_DEPLETION, "equals": True},
        },
        CONF_DEPLETION_INTEGRATOR: {
            **common,
            "order": 82,
            "showWhen": {"setting": CONF_DEPLETION, "equals": True},
        },
        CONF_DEPLETION_RATES: {
            **common,
            "order": 84,
            "showWhen": {"setting": CONF_DEPLETION, "equals": True},
        },
        CONF_DEPLETION_GROUPS: {
            **common,
            "order": 85,
            "showWhen": {"setting": CONF_DEPLETION_RATES, "equals": "flux"},
        },
        CONF_PYTHON: {
            **common,
            "order": 83,
            "widget": "path",
            "showWhen": {"setting": CONF_DEPLETION, "equals": True},
        },
    }
    return editors
