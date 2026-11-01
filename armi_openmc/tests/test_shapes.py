import math

from armi_openmc.shapes import Outline


def test_concentricCircles():
    assert Outline("circle", (1.0,)).contains(Outline("circle", (0.9,)))
    assert not Outline("circle", (0.9,)).contains(Outline("circle", (1.0,)))


def test_pinInsideFlatTopHexagonDependsOnDirection():
    duct = Outline("hexagon", (10.0,), hexOrientation="x")  # flats at top and bottom
    # Along x the flat-top hexagon reaches its corner, 10/sqrt(3) = 5.77 cm out.
    assert duct.contains(Outline("circle", (1.0,), x=5.0, y=0.0))
    # Along y it reaches only its flat, 5 cm out.
    assert not duct.contains(Outline("circle", (1.0,), x=0.0, y=5.0))


def test_areas():
    assert math.isclose(Outline("hexagon", (2.0,)).area, 2.0 * math.sqrt(3.0))
    assert math.isclose(Outline("rectangle", (2.0, 3.0)).area, 6.0)
