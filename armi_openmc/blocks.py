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
* *Spherical blocks* (a pebble's unit cell) are concentric spheres, each a component.
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
from armi.nucDirectory import nuclideBases
from armi.reactor import grids
from armi.reactor.blocks import HexBlock
from armi.reactor.components import DerivedShape, Sphere

from armi_openmc import particles as particleFill
from armi_openmc import shapes
from armi_openmc.materials import MaterialFactory, mixNumberDensities

try:  # the ARMI fork with spherical (pebble) unit cells
    from armi.reactor.blocks import SphericalBlock
except ImportError:  # pragma: no cover - older ARMI
    SphericalBlock = ()

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


def isSpherical(b):
    return bool(SphericalBlock) and isinstance(b, SphericalBlock)


def blockOutline(b, hexOrientation="x"):
    """The footprint of a block: its hexagon or rectangle of pitch size, or its outer sphere."""
    if isSpherical(b):
        return _rounded(shapes.Outline("sphere", (float(b.getPitch()),)))
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

    def __init__(self, materials: MaterialFactory, hexOrientation="x", options=None):
        self.materials = materials
        self.hexOrientation = hexOrientation
        self.surfaces = SurfaceCache()
        self.particleModel = getattr(options, "particleModel", "homogenized")
        self.particles = (
            particleFill.ParticleFiller(materials, options) if options is not None else None
        )
        self._zOrigin = 0.0
        self._inclusions = []

    def build(self, b, heterogeneous=True, zOrigin=0.0) -> BlockModel:
        """
        The universe of block ``b``.

        Parameters
        ----------
        zOrigin : float
            Elevation of the block's bottom in the model, where explicit particles start.
        """
        self._zOrigin = zOrigin
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
        if hasattr(self.materials, "forBlock"):  # multigroup: macroscopic data of the whole block
            mat = self.materials.forBlock(b)
            fissionableMat = self.materials.isFissionable(mat)
        else:
            mat = self.materials.make(
                f"{b.getName()}", b.getNumberDensities(), b.getAverageTempInC()
            )
            fissionableMat = _isFissionable(mat)
        cell = openmc.Cell(name=f"{b.getName()}", fill=mat)
        universe = openmc.Universe(name=b.getName(), cells=[cell])
        fissionable = [cell] if fissionableMat else []
        return BlockModel(b, universe, False, fissionable)

    def _heterogeneous(self, b):
        if isSpherical(b):
            return self._spherical(b)
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
            cells.append(self._componentCell(b, c, region, fissionable))
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

        # Off-center components (a control channel beside the pins) are carved out of the interior.
        inclusionOutlines = []
        for c in self._inclusions:
            x, y, _ = c.spatialLocator.getLocalCoordinates()
            outer, inner = shapes.outlines(c, orientation, float(x), float(y))
            if interior is not None and not interior.contains(outer):
                raise Unrepresentable(f"{c.name} reaches outside the innermost shell")
            inclusionOutlines.append(outer)
            carve = self.surfaces.outside(outer)
            interiorRegion = interiorRegion & carve if interiorRegion is not None else carve
            region = self.surfaces.inside(outer)
            if inner is not None:
                region &= self.surfaces.outside(inner)
                cells.append(
                    openmc.Cell(
                        name=f"{b.getName()} {c.name} center",
                        fill=bgMaterial,
                        region=self.surfaces.inside(inner),
                    )
                )
            cells.append(self._componentCell(b, c, region, fissionable))

        model = BlockModel(b, None, True)
        if pinsAt:
            lattice = self._pinLattice(
                b, pinsAt, bgMaterial, interior, fissionable, model, inclusionOutlines
            )
            cells.append(
                openmc.Cell(name=f"{b.getName()} pins", fill=lattice, region=interiorRegion)
            )
        else:
            cells.append(
                openmc.Cell(name=f"{b.getName()} interior", fill=bgMaterial, region=interiorRegion)
            )
        return self._finish(b, cells, fissionable, model)

    def _spherical(self, b):
        """Concentric spheres, innermost first, each filling the space up to the next."""
        spheres = sorted(
            (c for c in b.iterComponents() if not particleFill.isParticleLayer(c)),
            key=lambda c: c.getDimension("od"),
        )
        cells, fissionable = [], []
        inner = None
        for c in spheres:
            if not isinstance(c, Sphere) or c.getDimension("mult") != 1:
                raise Unrepresentable(f"{c.name} is not a single sphere")
            iD, od = (c.getDimension("id") or 0.0), c.getDimension("od")
            previous = inner.size[0] if inner is not None else 0.0
            if abs(iD - previous) > 1e-6:
                raise Unrepresentable(f"{c.name} starts at {iD} cm, not at {previous} cm")
            outer = _rounded(shapes.Outline("sphere", (float(od),)))
            region = self.surfaces.inside(outer)
            if inner is not None:
                region &= self.surfaces.outside(inner)
            cells.append(self._componentCell(b, c, region, fissionable))
            inner = outer
        return self._finish(b, cells, fissionable)

    def _finish(self, b, cells, fissionable, model=None):
        model = model or BlockModel(b, None, True)
        model.universe = openmc.Universe(name=b.getName(), cells=cells)
        model.fissionableCells = fissionable
        return model

    # ----------------------------------------------------------------------------------------
    def _classify(self, b):
        """Split components into shells, pins (by grid position) and background."""
        shells, background, inclusions = [], [], []
        pinsAt = defaultdict(list)
        for c in b.iterComponents():
            if particleFill.isParticleLayer(c):
                continue  # drawn inside its host
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
            elif c.getDimension("mult") == 1 and _offCenter(c):
                inclusions.append(c)
            elif c.getDimension("mult") == 1:
                shells.append(c)
            else:
                runLog.warning(
                    f"OpenMC: {c} has multiplicity {c.getDimension('mult')} but no grid positions; "
                    "smearing it into the block background.",
                    single=True,
                )
                background.append(c)
        self._inclusions = inclusions
        return shells, pinsAt, background

    def _backgroundMaterial(self, b, background):
        parts = [(c.getNumberDensities(), c.getArea(), c.temperatureInC) for c in background]
        parts = [p for p in parts if p[1] > 0.0]
        if not parts:
            return None
        ndens, temp = mixNumberDensities(parts)
        names = "+".join(c.name for c in background)
        return self.materials.make(f"{b.getName()} {names}", ndens, temp)

    def _componentCell(self, b, c, region, fissionable):
        """The cell of one component; a particle host is filled with them, or their mixture."""
        layers = particleFill.particleLayers(c)
        if not layers:
            mat = self._material(b, c)
            cell = openmc.Cell(name=f"{b.getName()} {c.name}", fill=mat, region=region)
            if _isFissionable(mat):
                fissionable.append(cell)
            return cell
        why = (
            self.particles.whyNotExplicit(c)
            if self.particles is not None
            else "no particle options"
        )
        if self.particleModel == "explicit" and why is None:
            cell = openmc.Cell(
                name=f"{b.getName()} {c.name}", fill=self.particles.fill(b, c), region=region
            )
            if not isinstance(c, Sphere):  # a prism's particles stack up from the block bottom
                cell.translation = (0.0, 0.0, self._zOrigin)
            if any(_hasActinides(layer) for layer in layers):
                fissionable.append(cell)
            return cell
        if self.particleModel == "explicit":
            runLog.warning(f"OpenMC: smearing the particles in {c.name}: {why}", single=True)
        if self.particles is not None:
            mat = self.particles.homogenizedMaterial(b, c)
        else:
            mat = self.materials.make(
                f"{b.getName()} {c.name} homogenized",
                *mixNumberDensities(
                    [(x.getNumberDensities(), x.getArea(), x.temperatureInC) for x in [c, *layers]]
                ),
            )
        cell = openmc.Cell(name=f"{b.getName()} {c.name}", fill=mat, region=region)
        if _isFissionable(mat):
            fissionable.append(cell)
        return cell

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
    def _pinLattice(self, b, pinsAt, bgMaterial, interior, fissionable, model, inclusions=()):
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
            placed = shapes.Outline(*_shifted(pinUniverses[key][1], x, y))
            if interior is not None and not interior.contains(placed):
                raise Unrepresentable("its pin lattice reaches outside the innermost shell")
            if any(placed.overlaps(other) for other in inclusions):
                raise Unrepresentable(f"a pin at {ij} overlaps an off-center component")
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
            cells.append(self._componentCell(b, c, region, fissionable))
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


def _offCenter(c):
    """True for a component whose location puts it off the center of its block."""
    loc = c.spatialLocator
    if not isinstance(loc, grids.CoordinateLocation):
        return False
    x, y, _ = loc.getLocalCoordinates()
    return abs(float(x)) > 1e-9 or abs(float(y)) > 1e-9


def _shifted(outline, x, y):
    return (outline.kind, outline.size, x, y, outline.hexOrientation)


def _setNested(nested, index, value):
    target = nested
    for i in index[:-1]:
        target = target[i]
    target[index[-1]] = value


def _hasActinides(component):
    """True for an ARMI component holding any actinide."""
    for name, density in component.getNumberDensities().items():
        base = nuclideBases.byName.get(name)
        if density > 0.0 and base is not None and base.z >= 90:
            return True
    return False


def _isFissionable(mat):
    if mat is None:
        return False
    return any(openmc.data.zam(n)[0] >= 90 for n in mat.get_nuclides())
