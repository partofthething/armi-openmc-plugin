"""
OpenMC as ARMI's global flux solver.

With ``neutronicsKernel: OpenMC`` the interface builds a whole-core model at every time node,
runs it and puts keff, block power, power density and multigroup flux on the reactor. Unlike a
diffusion code it needs no cross section library and no uniform axial mesh: each assembly keeps
its own block boundaries, so ``detailedAxialExpansion`` cores run as they are.
"""

import os

from armi import interfaces, runLog
from armi.physics.neutronics.globalFlux import globalFluxInterface as gfi
from armi.physics.neutronics.settings import CONF_NEUTRONICS_KERNEL

from armi_openmc import runner
from armi_openmc.core import MGXS_FILE, CoreModelBuilder
from armi_openmc.options import OpenMCOptions
from armi_openmc.results import CoreResults
from armi_openmc.settings import KERNEL_NAME

ORDER = interfaces.STACK_ORDER.FLUX


def describeInterfaces(cs):
    if cs[CONF_NEUTRONICS_KERNEL] == KERNEL_NAME:
        return (OpenMCGlobalFluxInterface, {})
    return None


class OpenMCFluxOptions(gfi.GlobalFluxOptions):
    """Global flux options plus the OpenMC modeling options."""

    def __init__(self, label=None):
        super().__init__(label)
        self.openmc = OpenMCOptions()
        self.power = None

    def fromUserSettings(self, cs):
        super().fromUserSettings(cs)
        self.openmc = OpenMCOptions.fromSettings(cs)
        self.power = cs["power"]


class OpenMCExecuter(gfi.GlobalFluxExecuter):
    """
    Builds, runs and reads one whole-core OpenMC model.

    It keeps its own run sequence rather than the framework's: OpenMC needs no geometry
    transformation, and its run directory is a plain folder per label under ``openmcRunDir``.
    """

    def run(self):
        opts = self.options
        workdir = os.path.join(opts.openmc.runDir, opts.label or "openmc")
        os.makedirs(workdir, exist_ok=True)
        coreModel = CoreModelBuilder(self.r.core, opts.openmc).build()
        modelPath = os.path.join(workdir, "model.xml")
        previous = _read(modelPath)
        if coreModel.multigroup is not None:
            coreModel.model.materials.cross_sections = os.path.abspath(
                os.path.join(workdir, MGXS_FILE)
            )
        coreModel.model.export_to_model_xml(modelPath)
        if coreModel.multigroup is not None:
            # the data are part of the model: a changed library is a changed model
            digestPath = os.path.join(workdir, MGXS_FILE + ".sha1")
            previous = previous if _read(digestPath) == coreModel.multigroup.digest else None
            coreModel.multigroup.export(os.path.join(workdir, MGXS_FILE))
            with open(digestPath, "w") as f:
                f.write(coreModel.multigroup.digest)
        statepoint = None
        if opts.openmc.reuseStatepoint and previous is not None and previous == _read(modelPath):
            statepoint = runner.lastStatepoint(workdir)
            if statepoint:
                runLog.important(
                    f"OpenMC: reusing {statepoint}, whose model is identical to this one"
                )
        if statepoint is None:
            if opts.openmc.plotGeometry:
                runner.run(opts.openmc, workdir, args=["-p"], logName="plot.out")
            statepoint = runner.run(opts.openmc, workdir)
        results = CoreResults(coreModel, statepoint)
        mapper = OpenMCResultMapper(results, opts.power, r=self.r, cs=opts.cs)
        if opts.applyResultsToReactor:
            mapper.apply(self.r)
        return mapper


def _read(path):
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        return f.read()


class OpenMCResultMapper(gfi.GlobalFluxResultMapper):
    def __init__(self, results: CoreResults, power, r=None, cs=None):
        super().__init__(r=r, cs=cs)
        self.results = results
        self.power = power

    def getKeff(self):
        return self.results.keff

    def apply(self, reactor):
        self.r = reactor
        self.results.apply(reactor, self.power)
        self._updateDerivedParams()
        self.success = True


class OpenMCGlobalFluxInterface(gfi.GlobalFluxInterfaceUsingExecuters):
    """Runs OpenMC on the whole core at every time node."""

    name = "openmcFlux"

    @staticmethod
    def getOptionsCls():
        return OpenMCFluxOptions

    @staticmethod
    def getExecuterCls():
        return OpenMCExecuter

    def interactBOL(self):
        super().interactBOL()
        opts = OpenMCOptions.fromSettings(self.cs)
        runLog.info(
            f"OpenMC global flux: {opts.particles} particles x {opts.batches} batches "
            f"({opts.inactive} inactive), {opts.blockModel} blocks, runner {opts.runner}"
        )
