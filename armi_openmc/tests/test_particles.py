"""Explicit and homogenized particle fuel in block universes (needs ARMI with particles)."""

import inspect
import math

import numpy as np
import openmc
import pytest

from armi_openmc.blocks import BlockBuilder
from armi_openmc.materials import MaterialFactory, NuclideMapper
from armi_openmc.options import OpenMCOptions
from armi_openmc.tests.conftest import crossSections, needsData

particles = pytest.importorskip("armi.reactor.components.particles")

ELEMENT = inspect.cleandoc(
    """
    nuclide flags:
        U235: {burn: false, xs: true}
        U238: {burn: false, xs: true}
        O: {burn: false, xs: true}
        C: {burn: false, xs: true}
        SI: {burn: false, xs: true}
        HE4: {burn: false, xs: true}
    particles:
        TRISO:
            kernel: {material: UCO, od: 0.0425, flags: fuel depletable}
            buffer: {material: Graphite, od: 0.0625}
            inner PyC: {material: Graphite, od: 0.0695}
            SiC: {material: SiC, od: 0.0765}
            outer PyC: {material: Graphite, od: 0.0845}
    blocks:
        fuel: &block_fuel
            grid name: element
            fuel compact:
                shape: Circle
                material: Graphite
                Tinput: 20.0
                Thot: 20.0
                od: 1.245
                latticeIDs: [F]
                particles: TRISO
                packingFraction: 0.35
            coolant:
                shape: Circle
                material: Helium
                Tinput: 20.0
                Thot: 20.0
                od: 1.588
                latticeIDs: [C]
            channel:
                shape: Circle
                material: Helium
                Tinput: 20.0
                Thot: 20.0
                od: 1.0
                origin: [3.0, 3.0]
            graphite:
                shape: DerivedShape
                material: Graphite
                Tinput: 20.0
                Thot: 20.0
            outline:
                shape: Hexagon
                material: Void
                Tinput: 20.0
                Thot: 20.0
                ip: 10.0
                op: 10.0
    assemblies:
        column:
            specifier: F
            blocks: [*block_fuel]
            height: [10.0]
            axial mesh points: [1]
            xs types: [A]
    grids:
        element:
            geom: hex_corners_up
            symmetry: full
            lattice pitch:
                hex: 1.8796
            grid contents:
                [0, 0]: C
                [1, 0]: F
                [0, 1]: F
                [-1, 1]: F
                [-1, 0]: F
                [0, -1]: F
                [1, -1]: F
    """
)


@pytest.fixture(scope="module")
def element():
    from armi import settings
    from armi.reactor import blueprints

    return blueprints.Blueprints.load(ELEMENT).constructAssem(settings.Settings(), "column")[0]


def _builder(model, tmp_path):
    opts = OpenMCOptions(crossSections=crossSections(), runDir=str(tmp_path), particleModel=model)
    return BlockBuilder(MaterialFactory(NuclideMapper(crossSections())), "x", options=opts)


@needsData
def test_explicitParticlesFillTheCompacts(element, tmp_path):
    model = _builder("explicit", tmp_path).build(element, zOrigin=0.0)
    assert model.heterogeneous
    trisos = [
        c for c in model.universe.get_all_cells().values() if isinstance(c, openmc.model.TRISO)
    ]
    host = element.getComponentByName("fuel compact")
    perCompact = host.particleLayers[0].getDimension("mult") / host.getDimension("mult") * 10.0
    # two 5 cm segments of the same packing: each TRISO cell is one particle of a segment, or a
    # copy where a particle straddles mesh elements
    assert len(trisos) >= round(perCompact / 2)
    # the compact cell is translated to the block's bottom and holds actinides
    compacts = [c for c in model.fissionableCells if "fuel compact" in c.name]
    assert compacts and all(np.allclose(c.translation, (0, 0, 0)) for c in compacts)
    # packings are cached for the next run
    assert list((tmp_path / "particle-packings").glob("*.npy"))


@needsData
def test_homogenizedConservesTheCompact(element, tmp_path):
    builder = _builder("homogenized", tmp_path)
    model = builder.build(element)
    host = element.getComponentByName("fuel compact")
    mixed = next(m for m in builder.materials.materials if "homogenized" in m.name)
    outline = host.getArea() + sum(c.getArea() for c in host.particleLayers)
    kernel = host.particleLayers[0]
    expected = kernel.getNumberDensity("U235") * kernel.getArea() / outline
    got = mixed.get_nuclide_atom_densities()["U235"]
    assert got == pytest.approx(expected, rel=1e-9)
    assert not any(
        isinstance(c, openmc.model.TRISO) for c in model.universe.get_all_cells().values()
    )


@needsData
def test_offCenterChannelIsCarvedOut(element, tmp_path):
    model = _builder("homogenized", tmp_path).build(element)
    channel = [c for c in model.universe.cells.values() if c.name.endswith("channel")]
    assert len(channel) == 1
    # the channel's cylinder sits at its origin with its radius
    (surface,) = channel[0].region.get_surfaces().values()
    assert (surface.x0, surface.y0, surface.r) == pytest.approx((3.0, 3.0, 0.5))
    # and the pin lattice cell excludes it
    pins = next(c for c in model.universe.cells.values() if c.name.endswith(" pins"))
    assert surface.id in pins.region.get_surfaces()
    assert math.isclose(element.getComponentByName("channel").getArea(), math.pi / 4.0)


@pytest.fixture(scope="module")
def pebble():
    from armi import settings
    from armi.reactor import blueprints

    bpTests = pytest.importorskip("armi.reactor.blueprints.tests.test_particleBlueprints")
    if not hasattr(bpTests, "PEBBLE_BLUEPRINT"):
        pytest.skip("ARMI without spherical blocks")
    packing = bpTests.PEBBLE_PARTICLES * particles.sphereVolume(0.09154) / particles.sphereVolume(5)
    text = bpTests.PEBBLE_BLUEPRINT.replace("PACKING", repr(packing))
    return blueprints.Blueprints.load(text).constructAssem(settings.Settings(), "pebble")[0]


@needsData
def test_pebbleCellIsConcentricSpheres(pebble, tmp_path):
    model = _builder("homogenized", tmp_path).build(pebble)
    radii = sorted(
        s.r for c in model.universe.cells.values() for s in c.region.get_surfaces().values()
    )
    assert sorted(set(radii)) == pytest.approx([2.5, 3.0, 3.7959])
    assert [c.name.split()[-1] for c in model.fissionableCells] == ["zone"]


@needsData
def test_explicitParticlesFillThePebble(pebble, tmp_path):
    model = _builder("explicit", tmp_path).build(pebble)
    (packing,) = (tmp_path / "particle-packings").glob("*.npy")
    centers = np.load(packing)
    assert len(centers) == 9394  # 9,393.6 particles, rounded
    assert np.linalg.norm(centers, axis=1).max() <= 2.5 - 0.09154 / 2 + 1e-9
    zone = next(c for c in model.fissionableCells if c.name.endswith("fuel zone"))
    assert zone.translation is None


@needsData
def test_pebbleCellHasAWhiteBoundary(pebble, tmp_path):
    from armi_openmc.mgxs import XSModelBuilder

    opts = OpenMCOptions(crossSections=crossSections(), runDir=str(tmp_path))
    xs = XSModelBuilder(opts).build(pebble, "AA", ["U235", "U238", "C"])
    surfaces = xs.model.geometry.get_all_surfaces().values()
    boundaries = {
        (round(s.r, 6), s.boundary_type) for s in surfaces if isinstance(s, openmc.Sphere)
    }
    assert [r for r, kind in boundaries if kind == "white"] == [3.7959]
