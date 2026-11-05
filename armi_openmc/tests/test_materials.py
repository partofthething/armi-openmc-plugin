import pytest

from armi_openmc.materials import MaterialFactory, NuclideMapper, mixNumberDensities
from armi_openmc.tests.conftest import crossSections, needsData


@pytest.fixture(scope="module")
def mapper():
    return NuclideMapper(crossSections())


@needsData
def test_isotopesAndElements(mapper):
    assert mapper.expand("U235") == [("U235", 1.0)]
    assert mapper.expand("AM242M") == [("Am242_m1", 1.0)]
    iron = dict(mapper.expand("FE"))
    assert set(iron) == {"Fe54", "Fe56", "Fe57", "Fe58"}
    assert sum(iron.values()) == pytest.approx(1.0)
    assert mapper.expand("LFP35") == []


@needsData
def test_waterGetsThermalScattering(mapper):
    mat = MaterialFactory(mapper).make("water", {"H1": 0.0667, "O16": 0.0334}, 20.0)
    assert mat._sab and mat._sab[0][0] == "c_H_in_H2O"
    assert mat.temperature == pytest.approx(293.15)


@needsData
def test_voidIsNone(mapper):
    assert MaterialFactory(mapper).make("gap", {}, 20.0) is None


def test_naturalEvaluationCarriesMissingIsotopes(tmp_path):
    """A library with only natural carbon (ENDF/B-VII.1's C0) still gets the carbon."""
    library = tmp_path / "cross_sections.xml"
    library.write_text(
        '<?xml version="1.0"?>\n<cross_sections>\n'
        '  <library materials="C0" path="C0.h5" type="neutron" />\n'
        '  <library materials="Fe56" path="Fe56.h5" type="neutron" />\n'
        "</cross_sections>\n"
    )
    mapper = NuclideMapper(str(library))
    assert mapper.expand("C12") == [("C0", 1.0)]
    assert mapper.expand("C") == [("C0", 1.0)]
    assert mapper.expand("FE56") == [("Fe56", 1.0)]
    assert mapper.expand("U235") == []


def test_mixByVolume():
    mixed, temp = mixNumberDensities([({"A": 1.0}, 1.0, 100.0), ({"A": 0.0, "B": 2.0}, 3.0, 200.0)])
    assert mixed == {"A": 0.25, "B": 1.5}
    assert temp == pytest.approx(175.0)
