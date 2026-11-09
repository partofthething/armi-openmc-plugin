"""Whole-core models of the example reactors (built, not run)."""

import openmc
from armi.reactor.flags import Flags

from armi_openmc.core import CoreModelBuilder
from armi_openmc.tests.conftest import needsData


@needsData
def test_c5g7QuarterStaysAQuarter(c5g7, options):
    model = CoreModelBuilder(c5g7.core, options).build()
    assert model.cut
    assert model.domainFraction == 0.25
    assert len(model.blockCells) == 36
    assert set(model.copies.values()) == {1}
    root = next(iter(model.model.geometry.root_universe.cells.values()))
    assert isinstance(root.fill, openmc.RectLattice)
    assert root.fill.shape == (3, 3)


@needsData
def test_c5g7CanUnfold(c5g7, options):
    options.symmetry = "full"
    model = CoreModelBuilder(c5g7.core, options).build()
    assert not model.cut
    assert max(model.copies.values()) == 4


@needsData
def test_c5g7Slice(c5g7, options):
    options.axialSlice = 30.0
    model = CoreModelBuilder(c5g7.core, options).build()
    assert len(model.blockCells) == 9


@needsData
def test_fftfFullHexCore(fftf, options):
    options.heterogeneousFlags = ["fuel"]
    model = CoreModelBuilder(fftf.core, options).build()
    assert not model.cut
    root = next(iter(model.model.geometry.root_universe.cells.values()))
    assert isinstance(root.fill, openmc.HexLattice)
    assert root.fill.orientation == "y"
    het = [b for b, m in model.blockModels.items() if m.heterogeneous and m.pinLattice is not None]
    assert het and all(b.hasFlags(Flags.FUEL) for b in het)


@needsData
def test_thirdCoreHexUnfolds(options):
    from armi.testing import loadTestReactor

    _o, r = loadTestReactor()  # the third-core sodium hex reactor
    options.blockModel = "homogenized"
    model = CoreModelBuilder(r.core, options).build()
    assert not model.cut
    copies = {b.parent.spatialLocator.getRingPos(): n for b, n in model.copies.items()}
    assert copies[(1, 1)] == 1
    assert all(n == 3 for (ring, _), n in copies.items() if ring > 1)


@needsData
def test_emptyPositionsHoldTheOuterFill(fftf, options):
    """Lattice positions without an assembly are what lies outside the assemblies, not void."""
    options.blockModel = "homogenized"
    options.outerFill = "Outer Radial Shield"
    model = CoreModelBuilder(fftf.core, options).build()
    root = next(iter(model.model.geometry.root_universe.cells.values()))
    universes = [u for ring in root.fill.universes for u in ring]
    assert not any(u.name == "void" for u in universes)
    assert any(u.name == "outer fill" for u in universes)
    assert root.fill.outer.name == "outer fill"


@needsData
def test_cartesianCoreInACylinder(c5g7, options):
    """A square lattice inside a larger cylinder, the space between filled (a pool, a drum)."""
    import pytest

    options.symmetry = "full"
    options.blockModel = "homogenized"
    options.coreRadius = 200.0
    fill = next(iter(c5g7.core.iterBlocks())).getType()
    options.outerFill = fill
    model = CoreModelBuilder(c5g7.core, options).build()
    root = next(iter(model.model.geometry.root_universe.cells.values()))
    cylinders = [s for s in root.region.get_surfaces().values() if isinstance(s, openmc.ZCylinder)]
    assert [c.r for c in cylinders] == [200.0]
    assert root.fill.outer.name == "outer fill"
    options.coreRadius = 10.0
    with pytest.raises(ValueError, match="cuts through"):
        CoreModelBuilder(c5g7.core, options).build()
