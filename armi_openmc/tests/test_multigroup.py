"""Tests for the multigroup mode: macroscopic OpenMC data from an ARMI ISOTXS library."""

import types

import numpy as np
import pytest
from scipy import sparse

from armi_openmc import multigroup


def _nuclide(total, fission, nu, chi, capture, n2n, scatter0, scatter1, efiss=0.0):
    """A stand-in for an ARMI library nuclide (micros as an ISOTXS gives them)."""
    micros = types.SimpleNamespace(
        total=np.asarray(total, float).reshape(-1, 1),
        fission=np.asarray(fission, float),
        neutronsPerFission=np.asarray(nu, float),
        chi=None if chi is None else np.asarray(chi, float),
        nGamma=np.asarray(capture, float),
        nalph=None,
        np=None,
        nd=None,
        nt=None,
        n2n=np.asarray(n2n, float),
        totalScatter=sparse.csc_matrix(np.asarray(scatter0, float)),
        higherOrderScatter={1: sparse.csc_matrix(np.asarray(scatter1, float))},
    )
    return types.SimpleNamespace(micros=micros, isotxsMetadata={"efiss": efiss})


class FakeLibrary(dict):
    numGroups = 2
    neutronEnergyUpperBounds = np.array([1.0e7, 1.0e5])


class FakeBlock:
    def __init__(self, name, densities, suffix="AA"):
        self.name, self.densities, self.suffix = name, densities, suffix

    def getName(self):
        return self.name

    def getMicroSuffix(self):
        return self.suffix

    def getNumberDensities(self):
        return self.densities


@pytest.fixture
def library():
    lib = FakeLibrary()
    # [to, from]: group 1 scatters 2.0 to itself, 1.0 down, plus 0.2 of (n,2n) counted twice
    lib["U235AA"] = _nuclide(
        total=[4.0, 6.0],
        fission=[1.0, 3.0],
        nu=[2.5, 2.4],
        chi=[0.9, 0.1],
        capture=[0.2, 0.5],
        n2n=[0.2, 0.0],
        scatter0=[[2.2, 0.0], [1.2, 2.5]],
        scatter1=[[0.3, 0.0], [0.1, 0.2]],
    )
    lib["NA23AA"] = _nuclide(
        total=[3.0, 4.0],
        fission=[0.0, 0.0],
        nu=[0.0, 0.0],
        chi=None,
        capture=[0.01, 0.02],
        n2n=[0.0, 0.0],
        scatter0=[[2.5, 0.0], [0.49, 3.98]],
        scatter1=[[0.2, 0.0], [0.0, 0.1]],
    )
    return lib


def test_macroscopicData(library):
    factory = multigroup.MultigroupMaterialFactory(library)
    mat = factory.forBlock(FakeBlock("fuel", {"U235": 0.01, "NA23": 0.02, "XX99": 1e-3}))
    assert factory.isFissionable(mat)
    (xs,) = factory.library.xsdatas
    T = multigroup.DATA_TEMPERATURE
    total = np.asarray(xs.total[0])
    absorption = np.asarray(xs.absorption[0])
    np.testing.assert_allclose(total, 0.01 * np.array([4.0, 6.0]) + 0.02 * np.array([3.0, 4.0]))
    np.testing.assert_allclose(
        absorption, 0.01 * np.array([1.2, 3.5]) + 0.02 * np.array([0.01, 0.02])
    )
    scatter = np.asarray(xs.scatter_matrix[0])  # [from, to, order]
    assert scatter.shape == (2, 2, 2)
    np.testing.assert_allclose(
        scatter[0, :, 0], 0.01 * np.array([2.2, 1.2]) + 0.02 * np.array([2.5, 0.49])
    )
    np.testing.assert_allclose(
        scatter[0, :, 1], 0.01 * np.array([0.3, 0.1]) + 0.02 * np.array([0.2, 0.0])
    )
    # multiplicity: production / (total - absorption), so collisions balance
    production = scatter[:, :, 0].sum(axis=1)
    multiplicity = np.asarray(xs.multiplicity_matrix[0])[:, 0]
    np.testing.assert_allclose(production / multiplicity, total - absorption)
    # nu-fission keeps each nuclide's spectrum: [from, to] = N nu sigma_f chi
    nuFission = np.asarray(xs.nu_fission[0])
    np.testing.assert_allclose(
        nuFission, 0.01 * np.outer([2.5, 2.4] * np.array([1.0, 3.0]), [0.9, 0.1])
    )
    np.testing.assert_allclose(
        np.asarray(xs.kappa_fission[0]), 0.01 * np.array([1.0, 3.0]) * multigroup.DEFAULT_KAPPA_EV
    )
    assert T == xs.temperatures[0]
    assert "XX99AA" in factory._missing


def test_identicalBlocksShareData(library):
    factory = multigroup.MultigroupMaterialFactory(library)
    a = factory.forBlock(FakeBlock("a", {"NA23": 0.02}))
    b = factory.forBlock(FakeBlock("b", {"NA23": 0.02}))
    c = factory.forBlock(FakeBlock("c", {"NA23": 0.03}))
    assert len(factory.library.xsdatas) == 2
    assert a is not b and a._macroscopic == b._macroscopic != c._macroscopic
    assert not factory.isFissionable(a)


def test_energyBoundsAscending(library):
    bounds = multigroup.energyBounds(library)
    np.testing.assert_allclose(bounds, [multigroup.LOWEST_ENERGY_EV, 1.0e5, 1.0e7])


def test_digestFollowsTheData(library):
    one = multigroup.MultigroupMaterialFactory(library)
    one.forBlock(FakeBlock("a", {"NA23": 0.02}))
    two = multigroup.MultigroupMaterialFactory(library)
    two.forBlock(FakeBlock("a", {"NA23": 0.02}))
    three = multigroup.MultigroupMaterialFactory(library)
    three.forBlock(FakeBlock("a", {"NA23": 0.021}))
    assert one.digest == two.digest != three.digest


def test_exportsAnOpenMCLibrary(library, tmp_path):
    import h5py

    factory = multigroup.MultigroupMaterialFactory(library)
    factory.forBlock(FakeBlock("fuel", {"U235": 0.01, "NA23": 0.02}))
    path = tmp_path / "mgxs.h5"
    factory.export(path)
    with h5py.File(path) as f:
        assert f.attrs["energy_groups"] == 2
        assert list(f.keys()) == ["AA-1"]


def test_isoaaFixture():
    """The real ISOTXS fixture of ARMI's tests, on its test fuel block."""
    from armi.nuclearDataIO.cccc import isotxs
    from armi.testing.singleBlocks import buildSimpleFuelHexBlock
    from armi.tests import ISOAA_PATH

    lib = isotxs.readBinary(ISOAA_PATH)
    block = buildSimpleFuelHexBlock()
    block.p.xsType = "A"
    factory = multigroup.MultigroupMaterialFactory(lib)
    mat = factory.forBlock(block)
    (xs,) = factory.library.xsdatas
    total = np.asarray(xs.total[0])
    assert total.shape == (lib.numGroups,) and (total > 0).all()
    assert (np.asarray(xs.absorption[0]) <= total).all()
    assert factory.isFissionable(mat)


def test_emptyGroupsGetDistinctBounds():
    lib = FakeLibrary()
    lib.neutronEnergyUpperBounds = np.array([1.0e7, 1.0e5, 1.0e5, 1.0e3])
    bounds = multigroup.energyBounds(lib)
    assert len(bounds) == 5
    assert (np.diff(bounds) > 0).all()
