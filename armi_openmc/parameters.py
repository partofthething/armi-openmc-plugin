"""Reactor parameters the OpenMC plugin adds."""

from armi.reactor import parameters
from armi.reactor.blocks import Block
from armi.reactor.parameters import ParamLocation
from armi.reactor.reactors import Core
from armi.utils import units


def getParameterDefinitions():
    corePDefs = parameters.ParameterDefinitionCollection()
    with corePDefs.createBuilder(
        default=0.0, location=ParamLocation.AVERAGE, categories=["neutronics"]
    ) as pb:
        pb.defParam(
            "openmcKeffStdDev",
            units=units.UNITLESS,
            description="One-sigma statistical uncertainty of the OpenMC eigenvalue",
        )

    blockPDefs = parameters.ParameterDefinitionCollection()
    with blockPDefs.createBuilder(
        default=None, location=ParamLocation.AVERAGE, categories=["neutronics"]
    ) as pb:
        pb.defParam(
            "openmcPinPowers",
            units=units.WATTS,
            description="Fission power of each fissionable pin, ordered like openmcPinIndices",
            saveToDB=True,
        )
        pb.defParam(
            "openmcPinIndices",
            units=units.UNITLESS,
            description="Block grid indices (i, j) of the pins in openmcPinPowers",
            saveToDB=True,
        )

    return {Core: corePDefs, Block: blockPDefs}
