"""
Surfaces for ARMI component shapes.

An ARMI component with a closed 2D shape (circle, hexagon, square, rectangle) is a prism along z
with an outer boundary and, when hollow, an inner one. :py:class:`Outline` captures one of those
boundaries so that it can become an OpenMC region and be compared with others for nesting.
"""

import math
from dataclasses import dataclass

import openmc
import openmc.model
from armi.reactor.components import basicShapes

SQRT3 = math.sqrt(3.0)

#: ARMI shape class -> (outer dimension names, inner dimension names)
_DIMS = {
    basicShapes.Circle: (("od",), ("id",)),
    basicShapes.Hexagon: (("op",), ("ip",)),
    basicShapes.Square: (("widthOuter",), ("widthInner",)),
    basicShapes.Rectangle: (("lengthOuter", "widthOuter"), ("lengthInner", "widthInner")),
    basicShapes.SolidRectangle: (("lengthOuter", "widthOuter"), ()),
}


def isSolidShape(component):
    """True for a component whose cross section OpenMC can draw exactly."""
    return type(component) in _DIMS


@dataclass(frozen=True)
class Outline:
    """One closed 2D boundary: a circle, a hexagon, or a rectangle, centered at (x, y)."""

    kind: str  # "circle", "hexagon" or "rectangle"
    size: tuple  # diameter | flat-to-flat | (x width, y width)
    x: float = 0.0
    y: float = 0.0
    hexOrientation: str = "x"  # OpenMC HexagonalPrism orientation; 'x' = flats on top

    @property
    def inscribedRadius(self):
        if self.kind == "rectangle":
            return min(self.size) / 2.0
        return self.size[0] / 2.0

    @property
    def circumscribedRadius(self):
        if self.kind == "circle":
            return self.size[0] / 2.0
        if self.kind == "hexagon":
            return self.size[0] / SQRT3
        return math.hypot(*self.size) / 2.0

    @property
    def area(self):
        if self.kind == "circle":
            return math.pi * self.size[0] ** 2 / 4.0
        if self.kind == "hexagon":
            return SQRT3 / 2.0 * self.size[0] ** 2
        return self.size[0] * self.size[1]

    def contains(self, other, tol=1e-6):
        """
        True if ``other`` lies entirely inside this outline (touching allowed).

        Concentric outlines of one kind compare dimensions. Otherwise ``other`` is replaced by
        its circumscribed circle, which is exact for round pins and conservative for the rest,
        and that circle is tested against each flat of this outline.
        """
        if (self.x, self.y) == (other.x, other.y) and self.kind == other.kind:
            if self.kind != "hexagon" or self.hexOrientation == other.hexOrientation:
                return all(o <= s + tol for o, s in zip(other.size, self.size, strict=True))
        dx, dy = other.x - self.x, other.y - self.y
        r = other.circumscribedRadius
        if self.kind == "circle":
            return math.hypot(dx, dy) + r <= self.size[0] / 2.0 + tol
        if self.kind == "rectangle":
            return (
                abs(dx) + r <= self.size[0] / 2.0 + tol and abs(dy) + r <= self.size[1] / 2.0 + tol
            )
        apothem = self.size[0] / 2.0
        first = 30.0 if self.hexOrientation == "x" else 0.0  # flat normals of the hexagon
        for k in range(6):
            angle = math.radians(first + 60.0 * k)
            if dx * math.cos(angle) + dy * math.sin(angle) + r > apothem + tol:
                return False
        return True

    def overlaps(self, other, tol=1e-6):
        """True if the two outlines' circumscribed circles intersect (conservative for polygons)."""
        distance = math.hypot(other.x - self.x, other.y - self.y)
        return distance < self.circumscribedRadius + other.circumscribedRadius - tol

    def same(self, other, tol=1e-6):
        return self.contains(other, tol) and other.contains(self, tol)

    def interior(self, boundaryType="transmission"):
        """The OpenMC region inside this outline (an infinite prism along z)."""
        if self.kind == "circle":
            surf = openmc.ZCylinder(
                x0=self.x, y0=self.y, r=self.size[0] / 2.0, boundary_type=boundaryType
            )
            return -surf
        if self.kind == "hexagon":
            prism = openmc.model.HexagonalPrism(
                edge_length=self.size[0] / SQRT3,
                orientation=self.hexOrientation,
                origin=(self.x, self.y),
                boundary_type=boundaryType,
            )
            return -prism
        prism = openmc.model.RectangularPrism(
            width=self.size[0],
            height=self.size[1],
            axis="z",
            origin=(self.x, self.y),
            boundary_type=boundaryType,
        )
        return -prism


def _kind(component):
    if isinstance(component, basicShapes.Circle):
        return "circle"
    if isinstance(component, basicShapes.Hexagon):
        return "hexagon"
    return "rectangle"


def outlines(component, hexOrientation="x", x=0.0, y=0.0):
    """
    The (outer, inner) outlines of a component at (x, y); inner is None for a solid shape.

    Dimensions are the hot ones ARMI uses for the component's area, so the OpenMC volumes
    match the ARMI number densities.
    """
    outerNames, innerNames = _DIMS[type(component)]
    kind = _kind(component)

    def dims(names):
        values = tuple(float(component.getDimension(n)) for n in names)
        if kind == "rectangle" and len(values) == 1:
            values = (values[0], values[0])
        return values

    outer = Outline(kind, dims(outerNames), x, y, hexOrientation)
    innerSize = dims(innerNames) if innerNames else (0.0,)
    inner = Outline(kind, innerSize, x, y, hexOrientation) if min(innerSize) > 0.0 else None
    return outer, inner
