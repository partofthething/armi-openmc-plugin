"""
OpenMC universes for ARMI blocks.

A block becomes a universe that is infinite in z and in the plane; whoever places it (a core
lattice, an infinite-lattice cell) clips it to the block's footprint. Two representations:

**Homogenized.** One cell filled with the block's volume-averaged number densities.

**Heterogeneous.** The components are laid out as ARMI describes them:

* *Pins* are components on the block's spatial grid (a ``MultiIndexLocation``). The components at
  each grid position, sorted by size, make concentric layers of a pin universe; the pin universes
  fill a rectangular or hexagonal lattice with the block grid's pitch and orientation.
* *Shells* are single centered shapes (ducts, cans, liners, an intercoolant hexagon), nested by
  size around the pin lattice.
* Everything else fills the space that is left. ``DerivedShape`` coolant is what ARMI already
  defines that way. Components OpenMC cannot draw exactly, such as helical wire wraps, are smeared
  into that same background by area, which is what they physically displace.

Every layer uses the hot dimensions ARMI computes its areas from, so the volume of each OpenMC
cell matches the ARMI component and the number densities carry over unchanged. A block whose
components cannot be laid out consistently (overlapping layers, pins that do not fit the pitch, a
pin grid with no room inside the innermost shell) falls back to the homogenized representation
with a warning that says why.
"""

from collections import defaultdict
from dataclasses import dataclass, field

import openmc
import openmc.data
from armi import runLog
from armi.reactor import grids
from armi.reactor.blocks import HexBlock
from armi.reactor.components import DerivedShape

from armi_openmc import shapes
from armi_openmc.materials import MaterialFactory, mixNumberDensities

ROUND = 9  # decimals kept on dimensions so that equal ARMI dimensions give identical surfaces


class Unrepresentable(Exception):
    """A block whose components cannot be laid out as heterogeneous geometry."""


@dataclass
class BlockModel:
    """The OpenMC representation of one ARMI block."""

    block: object
    universe: openmc.Universe
    heterogeneous: bool
    #: Material cells holding fissionable material, for pin-level tallies.
    fissionableCells: list = field(default_factory=list)
    #: Pin lattice of the block, if it has one.
    pinLattice: openmc.Lattice | None = None
    #: {pin lattice universe index: ARMI grid indices (i, j)} of the pin lattice.
    pinIndices: dict = field(default_factory=dict)


class SurfaceCache:
    """
    One OpenMC surface per distinct outline.

    Neighboring cells must share the surface between them. Two coincident but separate surfaces
    leave a particle on the wrong side of one of them and lose it.
    """

    def __init__(self):
        self._surfaces = {}

    def get(self, outline: shapes.Outline):
        key = _rounded(outline)
        if key not in self._surfaces:
            self._surfaces[key] = key.interior()
        return self._surfaces[key]

    def inside(self, outline):
        return self.get(outline)

    def outside(self, outline):
        return ~self.get(outline)


def _rounded(outline):
    return shapes.Outline(
        outline.kind,
        tuple(round(s, ROUND) for s in outline.size),
        round(outline.x, ROUND),
        round(outline.y, ROUND),
        outline.hexOrientation,
    )


def hexOrientationOf(grid):
    """
    OpenMC orientation of the hexagons that tile a hex grid.

    ARMI's default hex grid (``cornersUp=False``) has flat-topped cells, neighbors straight above
    one another. OpenMC calls a flat-topped hexagon ``'x'`` (two sides parallel to x) and a lattice
    of them ``'y'`` (columns along y).
    """
    return "y" if getattr(grid, "cornersUp", False) else "x"


def hexLatticeOrientationOf(grid):
    return "x" if getattr(grid, "cornersUp", False) else "y"


def blockOutline(b, hexOrientation="x"):
    """The footprint of a block: its hexagon or rectangle of pitch size."""
    if isinstance(b, HexBlock):
        return _rounded(
            shapes.Outline("hexagon", (float(b.getPitch()),), hexOrientation=hexOrientation)
        )
    pitch = b.getPitch()
    if isinstance(pitch, (tuple, list)):
        size = (float(pitch[0]), float(pitch[1]))
    else:
        size = (float(pitch), float(pitch))
    return _rounded(shapes.Outline("rectangle", size))


class BlockBuilder:
    """
    Builds block universes.

    Parameters
    ----------
    materials : MaterialFactory
        Makes (and records) every material the blocks use.
    hexOrientation : str
        OpenMC orientation of hexagonal blocks in the lattice that holds them ('x' = flats up).
    """

    def __init__(self, materials: MaterialFactory, hexOrientation="x"):
        self.materials = materials
        self.hexOrientation = hexOrientation
        self.surfaces = SurfaceCache()

    def build(self, b, heterogeneous=True) -> BlockModel:
        if heterogeneous:
            try:
                return self._heterogeneous(b)
            except Unrepresentable as why:
                runLog.warning(
                    f"OpenMC: homogenizing {b.getType()} ({b.getName()}): {why}",
                    single=True,
                    label=f"openmc-homog-{b.getType()}",
                )
        return self._homogenized(b)

    # ----------------------------------------------------------------------------------------
    def _homogenized(self, b):
        mat = self.materials.make(f"{b.getName()}", b.getNumberDensities(), b.getAverageTempInC())
        cell = openmc.Cell(name=f"{b.getName()}", fill=mat)
        universe = openmc.Universe(name=b.getName(), cells=[cell])
        fissionable = [cell] if _isFissionable(mat) else []
        return BlockModel(b, universe, False, fissionable)

    def _heterogeneous(self, b):
        shells, pinsAt, background = self._classify(b)
        if not shells and not pinsAt:
            raise Unrepresentable("it has no components with a drawable shape")
        bgMaterial = self._backgroundMaterial(b, background)
        orientation = self.hexOrientation
        footprint = blockOutline(b, orientation)

        layers = self._nest(
            [(c, *shapes.outlines(c, orientation)) for c in shells], f"shells of {b.getType()}"
        )
        cells, fissionable = [], []

        # Shells from the inside out, with background in any gap between them.
        for k, (c, outer, inner) in enumerate(layers):
            isOutermost = k == len(layers) - 1
            region = None
            if not (isOutermost and outer.same(footprint)):
                region = self.surfaces.inside(outer)
            if inner is not None:
                region = (
                    self.surfaces.outside(inner) & region
                    if region
                    else self.surfaces.outside(inner)
                )
            mat = self._material(b, c)
            cells.append(openmc.Cell(name=f"{b.getName()} {c.name}", fill=mat, region=region))
            if _isFissionable(mat):
                fissionable.append(cells[-1])
            if k > 0:
                prevOuter = layers[k - 1][1]
                if inner is not None and not inner.same(prevOuter):
                    gap = self.surfaces.inside(inner) & self.surfaces.outside(prevOuter)
                    cells.append(
                        openmc.Cell(name=f"{b.getName()} gap", fill=bgMaterial, region=gap)
                    )
        if layers:
            lastOuter = layers[-1][1]
            if not lastOuter.same(footprint):
                cells.append(
                    openmc.Cell(
                        name=f"{b.getName()} outside",
                        fill=bgMaterial,
                        region=self.surfaces.outside(lastOuter),
                    )
                )

        # The space inside the innermost shell holds the pins, in background.
        interior = None
        if layers:
            interior = layers[0][2]
            if interior is None:
                if pinsAt:
                    raise Unrepresentable("its pins sit inside a solid shell")
                return self._finish(b, cells, fissionable)
        interiorRegion = self.surfaces.inside(interior) if interior is not None else None

        model = BlockModel(b, None, True)
        if pinsAt:
            lattice = self._pinLattice(b, pinsAt, bgMaterial, interior, fissionable, model)
            cells.append(
                openmc.Cell(name=f"{b.getName()} pins", fill=lattice, region=interiorRegion)
            )
        else:
            cells.append(
                openmc.Cell(name=f"{b.getName()} interior", fill=bgMaterial, region=interiorRegion)
            )
        return self._finish(b, cells, fissionable, model)

    def _finish(self, b, cells, fissionable, model=None):
        model = model or BlockModel(b, None, True)
        model.universe = openmc.Universe(name=b.getName(), cells=cells)
        model.fissionableCells = fissionable
        return model

    # ----------------------------------------------------------------------------------------
    def _classify(self, b):
        """Split components into shells, pins (by grid position) and background."""
        shells, background = [], []
        pinsAt = defaultdict(list)
        for c in b.iterComponents():
            if isinstance(c, DerivedShape) or not shapes.isSolidShape(c):
                background.append(c)
                continue
            area = c.getArea()
            if area < -1e-10:
                raise Unrepresentable(f"component {c.name} has a negative area ({area:.4g} cm^2)")
            if area <= 1e-10:
                continue  # zero-area bookkeeping components such as a pitch-defining void
            loc = c.spatialLocator
            if isinstance(loc, grids.MultiIndexLocation) and b.spatialGrid is not None:
                for sub in loc:
                    i, j, _k = sub.getCompleteIndices()
                    pinsAt[(i, j)].append(c)
            elif c.getDimension("mult") == 1:
                shells.append(c)
            else:
                runLog.warning(
                    f"OpenMC: {c} has multiplicity {c.getDimension('mult')} but no grid positions; "
                    "smearing it into the block background.",
                    single=True,
                )
                background.append(c)
        return shells, pinsAt, background

    def _backgroundMaterial(self, b, background):
        parts = [(c.getNumberDensities(), c.getArea(), c.temperatureInC) for c in background]
        parts = [p for p in parts if p[1] > 0.0]
        if not parts:
            return None
        ndens, temp = mixNumberDensities(parts)
        names = "+".join(c.name for c in background)
        return self.materials.make(f"{b.getName()} {names}", ndens, temp)

    def _material(self, b, c):
        return self.materials.make(
            f"{b.getName()} {c.name}", c.getNumberDensities(), c.temperatureInC
        )

    @staticmethod
    def _nest(layers, what):
        """Sort (component, outer, inner) by size and check each fits inside the next."""
        layers = sorted(layers, key=lambda t: t[1].area)
        for (c1, outer1, _), (c2, _outer2, inner2) in zip(layers, layers[1:], strict=False):
            if inner2 is None or not inner2.contains(outer1):
                raise Unrepresentable(f"{what}: {c1.name} does not fit inside {c2.name}")
        return layers

    # ----------------------------------------------------------------------------------------
    def _pinLattice(self, b, pinsAt, bgMaterial, interior, fissionable, model):
        grid = b.spatialGrid
        bgCell = openmc.Cell(name=f"{b.getName()} pin background", fill=bgMaterial)
        bgUniverse = openmc.Universe(cells=[bgCell])

        pinUniverses = {}
        coords = {}
        for ij, comps in pinsAt.items():
            key = tuple(sorted(id(c) for c in comps))
            if key not in pinUniverses:
                pinUniverses[key] = self._pinUniverse(b, comps, bgMaterial, fissionable)
            x, y, _z = grid[ij[0], ij[1], 0].getLocalCoordinates()
            coords[ij] = (float(x), float(y), key)

        if isinstance(grid, grids.HexGrid):
            lattice = self._hexLattice(b, grid, coords, bgUniverse)
            cellOutline = shapes.Outline(
                "hexagon", (grid.pitch,), hexOrientation=hexOrientationOf(grid)
            )
        else:
            lattice = self._rectLattice(b, grid, coords, bgUniverse)
            cellOutline = shapes.Outline("rectangle", tuple(lattice.pitch))
        lattice.outer = bgUniverse

        universes = lattice.universes
        for ij, (x, y, key) in coords.items():
            idx = lattice.find_element((x, y, 0.0))[0]
            uIndex = lattice.get_universe_index(idx)
            _setNested(universes, uIndex, pinUniverses[key][0])
            model.pinIndices[uIndex] = ij
            if interior is not None:
                placed = shapes.Outline(*_shifted(pinUniverses[key][1], x, y))
                if not interior.contains(placed):
                    raise Unrepresentable("its pin lattice reaches outside the innermost shell")
        for _universe, outer in pinUniverses.values():
            if not cellOutline.contains(outer):
                raise Unrepresentable(f"its pins ({outer.size}) are larger than the pin pitch")
        lattice.universes = universes
        model.pinLattice = lattice
        return lattice

    def _pinUniverse(self, b, comps, bgMaterial, fissionable):
        """Concentric layers of one kind of pin, centered at the origin, in background."""
        layers = self._nest(
            [(c, *shapes.outlines(c, self._pinHexOrientation(b))) for c in comps], "pin"
        )
        cells = []
        for k, (c, outer, inner) in enumerate(layers):
            region = self.surfaces.inside(outer)
            if inner is not None:
                region &= self.surfaces.outside(inner)
            mat = self._material(b, c)
            cells.append(openmc.Cell(name=f"{b.getName()} {c.name}", fill=mat, region=region))
            if _isFissionable(mat):
                fissionable.append(cells[-1])
            if k == 0 and inner is not None:
                cells.append(
                    openmc.Cell(
                        name=f"{b.getName()} pin center",
                        fill=bgMaterial,
                        region=self.surfaces.inside(inner),
                    )
                )
            if k > 0 and inner is not None and not inner.same(layers[k - 1][1]):
                gap = self.surfaces.inside(inner) & self.surfaces.outside(layers[k - 1][1])
                cells.append(
                    openmc.Cell(name=f"{b.getName()} pin gap", fill=bgMaterial, region=gap)
                )
        outermost = layers[-1][1]
        cells.append(
            openmc.Cell(
                name=f"{b.getName()} pin coolant",
                fill=bgMaterial,
                region=self.surfaces.outside(outermost),
            )
        )
        return openmc.Universe(cells=cells), outermost

    def _pinHexOrientation(self, b):
        return hexOrientationOf(b.spatialGrid) if b.spatialGrid is not None else self.hexOrientation

    @staticmethod
    def _rectLattice(b, grid, coords, bgUniverse):
        px, py = (float(p) for p in grid.pitch)
        xs = [c[0] for c in coords.values()]
        ys = [c[1] for c in coords.values()]
        nx = int(round((max(xs) - min(xs)) / px)) + 1
        ny = int(round((max(ys) - min(ys)) / py)) + 1
        lattice = openmc.RectLattice(name=f"{b.getName()} pins")
        lattice.pitch = (px, py)
        lattice.lower_left = (min(xs) - px / 2.0, min(ys) - py / 2.0)
        lattice.universes = [[bgUniverse] * nx for _ in range(ny)]
        return lattice

    @staticmethod
    def _hexLattice(b, grid, coords, bgUniverse):
        rings = max(grid.getRingPos(ij)[0] for ij in coords)
        x0, y0, _ = grid[0, 0, 0].getLocalCoordinates()
        lattice = openmc.HexLattice(name=f"{b.getName()} pins")
        lattice.orientation = hexLatticeOrientationOf(grid)
        lattice.center = (float(x0), float(y0))
        lattice.pitch = (float(grid.pitch),)
        lattice.universes = [[bgUniverse] * max(1, 6 * (r - 1)) for r in range(rings, 0, -1)]
        return lattice


def _shifted(outline, x, y):
    return (outline.kind, outline.size, x, y, outline.hexOrientation)


def _setNested(nested, index, value):
    target = nested
    for i in index[:-1]:
        target = target[i]
    target[index[-1]] = value


def _isFissionable(mat):
    if mat is None:
        return False
    return any(openmc.data.zam(n)[0] >= 90 for n in mat.get_nuclides())
