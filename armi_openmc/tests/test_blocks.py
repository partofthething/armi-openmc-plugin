"""Block universes from the example reactors."""

import openmc

from armi_openmc.blocks import BlockBuilder
from armi_openmc.materials import MaterialFactory, NuclideMapper
from armi_openmc.tests.conftest import crossSections, needsData


def _builder():
    return BlockBuilder(MaterialFactory(NuclideMapper(crossSections())))


def _block(r, blockType):
    return next(b for b in r.core.iterBlocks() if b.getType() == blockType)


@needsData
def test_c5g7AssemblyIsPinByPin(c5g7):
    model = _builder().build(_block(c5g7, "mox"))
    assert model.heterogeneous
    assert isinstance(model.pinLattice, openmc.RectLattice)
    assert model.pinLattice.shape == (17, 17)
    assert len(model.pinIndices) == 17 * 17
    # mox low/medium/high fuel and the fission chamber's trace U-235
    assert len(model.fissionableCells) == 4


@needsData
def test_fftfFuelIsAHexBundleInADuct(fftf):
    model = _builder().build(_block(fftf, "Inner Fuel Pin"))
    assert model.heterogeneous
    assert isinstance(model.pinLattice, openmc.HexLattice)
    assert model.pinLattice.num_rings == 9
    assert len(model.pinIndices) == 217


@needsData
def test_fftfAbsorberHasTwoDucts(fftf):
    model = _builder().build(_block(fftf, "Absorber Control"))
    assert model.heterogeneous
    assert len(model.pinIndices) == 61
    names = {c.name.split(" ", 2)[-1] for c in model.universe.cells.values()}
    assert any("innerduct" in n for n in names)


@needsData
def test_homogenizedOnRequest(fftf):
    model = _builder().build(_block(fftf, "Inner Fuel Pin"), heterogeneous=False)
    assert not model.heterogeneous
    assert len(model.universe.cells) == 1
