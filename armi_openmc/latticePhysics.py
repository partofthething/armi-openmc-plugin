"""
OpenMC as ARMI's lattice physics code (``xsKernel: OpenMC``).

For every representative block of the cross section group manager it builds an OpenMC model
(:py:mod:`armi_openmc.mgxs`), runs it, reduces the tallies to microscopic group constants and
gathers them into one ISOTXS library that ARMI and any diffusion or transport plugin can use.
"""

import os
import shutil

from armi import interfaces, nuclearDataIO, runLog
from armi.nuclearDataIO.cccc import isotxs
from armi.physics import neutronics
from armi.physics.neutronics.const import CONF_CROSS_SECTION
from armi.physics.neutronics.latticePhysics.latticePhysicsInterface import LatticePhysicsInterface
from armi.physics.neutronics.settings import CONF_XS_KERNEL
from armi.reactor import grids

from armi_openmc import mgxs, runner
from armi_openmc.options import OpenMCOptions
from armi_openmc.settings import KERNEL_NAME

ORDER = interfaces.STACK_ORDER.CROSS_SECTIONS


def describeInterfaces(cs):
    if cs[CONF_XS_KERNEL] == KERNEL_NAME:
        return (OpenMCLatticePhysicsInterface, {})
    return None


def generateLibrary(r, cs, blocks, options: OpenMCOptions = None, suffix="", xsSettings=None):
    """
    Run OpenMC on representative blocks and return an ISOTXS library for them.

    Usable without the interface stack, e.g. from a script or a test.

    Parameters
    ----------
    r : Reactor
    cs : Settings
    blocks : list of Block
        Representative blocks, one per XS ID.
    """
    options = options or OpenMCOptions.fromSettings(cs)
    xsSettings = xsSettings if xsSettings is not None else cs[CONF_CROSS_SECTION]
    byID = {b.getMicroSuffix(): b for b in blocks}
    grid = r.core.spatialGrid if r is not None else None
    builder = mgxs.XSModelBuilder(options, grid if isinstance(grid, grids.HexGrid) else None)
    nuclides = list(r.blueprints.allNuclidesInProblem) if r is not None else []
    results = []
    for xsID, b in byID.items():
        driver = None
        driverID = getattr(xsSettings.get(xsID), "driverID", None) if xsSettings else None
        if options.xsColorset and driverID:
            driver = byID.get(driverID)
            if driver is None:
                runLog.warning(
                    f"OpenMC XS: driver {driverID} of {xsID} has no representative block"
                )
        blockNuclides = set(nuclides) | {n for n, d in b.getNumberDensities().items() if d > 0}
        xsModel = builder.build(b, xsID, blockNuclides, driver)
        workdir = mgxs.runDirFor(options, xsID, suffix)
        os.makedirs(workdir, exist_ok=True)
        xsModel.model.export_to_model_xml(os.path.join(workdir, "model.xml"))
        runLog.info(f"OpenMC XS {xsID}: {b} ({', '.join(xsModel.notes)})")
        statepoint = runner.run(options, workdir)
        constants, info = mgxs.reduce(xsModel, statepoint)
        if info["keff"] is not None:
            runLog.info(f"OpenMC XS {xsID}: k = {info['keff']}")
        results.append((xsModel, constants, info))
    return mgxs.buildLibrary(results, options.groupStructure, builder.mapper)


class OpenMCLatticePhysicsInterface(LatticePhysicsInterface):
    """Makes the ISOTXS library with OpenMC."""

    name = "openmcXS"

    def _getExecutablePath(self):
        return shutil.which(self.cs["openmcExecutable"]) or self.cs["openmcExecutable"]

    def _generateXsLibrary(
        self, baseList, forceSerial, xsLibrarySuffix, blockList, writers=None, purgeFP=True
    ):
        blocks = list(blockList or self.getRepresentativeBlocks())
        if baseList:
            blocks = [b for b in blocks if b.getMicroSuffix() in baseList]
        lib = generateLibrary(self.r, self.cs, blocks, suffix=xsLibrarySuffix)
        mgxs.writeLibrary(lib, neutronics.ISOTXS)
        return lib

    def readExistingXSLibraries(self, cycle, node):
        name = nuclearDataIO.getExpectedISOTXSFileName(cycle, node)
        path = name if os.path.exists(name) else neutronics.ISOTXS
        self.r.core.lib = isotxs.readBinary(path)

    def makeCycleXSFilesAsBaseFiles(self, cycle, node):
        self._copyLibraryFilesForCycle(
            cycle, {neutronics.ISOTXS: nuclearDataIO.getExpectedISOTXSFileName(cycle, node)}
        )

    def clearXS(self):
        pass

    def _getProcessesPerNode(self):
        return 1
