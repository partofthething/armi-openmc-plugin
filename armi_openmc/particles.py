"""
Explicit particle fuel: TRISO particles randomly packed in their host.

ARMI describes dispersed particles as a host component (a compact, whose material is the matrix) and
one ``ParticleLayer`` component per spherical layer, counted per centimeter of block height (see
``armi.reactor.components.particles`` in the ARMI fork). Here those become OpenMC geometry:

* one particle universe of concentric spheres, kernel to outer coating;
* particle centers randomly packed (``openmc.model.pack_spheres``) inside the host at ARMI's
  particle count;
* the particles sorted onto a Cartesian mesh (``openmc.model.create_triso_lattice``), so that a
  neutron only ever checks the handful of particles in its own mesh element rather than every
  particle in the compact, with the matrix as background.

Packing a whole block height in one go (10^5 particles and more) is slow and buys nothing, so the
host is packed in segments of about ``openmcParticleSegmentLength`` and the segment is stacked up
the block. Each segment is a closed packing domain, like the individual compacts stacked in a real
fuel hole: particles do not cross segment faces.

A spherical host (the fuel zone of a pebble) is packed whole, centered on the origin.

Identical hosts (same dimensions, particle count and compositions) share one packed universe, so a
core of fresh fuel blocks packs a compact once, and packings are cached on disk by their parameters.
"""

import hashlib
import math
import os

import numpy as np
import openmc
import openmc.model
from armi import runLog

from armi_openmc.materials import mixNumberDensities

try:  # the ARMI fork with particle fuel
    from armi.reactor.components.particles import ParticleLayer
except ImportError:  # pragma: no cover - older ARMI without particles
    ParticleLayer = ()

#: Particles per Cartesian mesh element the lattice pitch aims for.
PARTICLES_PER_ELEMENT = 8


def isParticleLayer(component):
    return bool(ParticleLayer) and isinstance(component, ParticleLayer)


def particleLayers(component):
    """The particle layers dispersed in a component (empty for an ordinary component)."""
    return list(getattr(component, "particleLayers", None) or [])


class ParticleFiller:
    """
    Builds the fill of particle-bearing hosts, explicit or homogenized.

    Parameters
    ----------
    materials : MaterialFactory
    options : OpenMCOptions
    """

    def __init__(self, materials, options):
        self.materials = materials
        self.options = options
        self._universes = {}
        self._particleUniverses = {}

    # ----------------------------------------------------------------------------------------
    def homogenizedMaterial(self, b, host):
        """Matrix and particles smeared into one material, by volume."""
        parts = [(host.getNumberDensities(), host.getArea(), host.temperatureInC)]
        parts += [
            (c.getNumberDensities(), c.getArea(), c.temperatureInC) for c in particleLayers(host)
        ]
        ndens, temp = mixNumberDensities(parts)
        return self.materials.make(f"{b.getName()} {host.name} homogenized", ndens, temp)

    def whyNotExplicit(self, host):
        """None if the host can be packed explicitly, else the reason it cannot."""
        if type(host).__name__ not in ("Circle", "Sphere"):
            return (
                f"particles in a {type(host).__name__} are only drawn in solid circles or spheres"
            )
        if (host.getDimension("id") or 0.0) > 0.0:
            return "particles in a hollow host cannot be packed"
        return None

    def _materials(self, host):
        layers = particleLayers(host)
        layerMaterials = [
            self.materials.make(f"{c.name}", c.getNumberDensities(), c.temperatureInC, reuse=True)
            for c in layers
        ]
        matrix = self.materials.make(
            f"{host.name} matrix", host.getNumberDensities(), host.temperatureInC, reuse=True
        )
        radii = tuple(round(c.getDimension("od") / 2.0, 9) for c in layers)
        return radii, layerMaterials, matrix

    def fill(self, b, host):
        """
        A universe of the host's particles, packed explicitly, centered on the host's axis.

        The universe spans z from 0 to the block height; the cell it fills must be translated to
        the block's bottom.
        """
        if type(host).__name__ == "Sphere":
            return self._fillSphere(host)
        layers = particleLayers(host)
        height = b.getHeight()
        segments = max(1, round(height / self.options.particleSegmentLength))
        segment = height / segments
        radius = host.getDimension("od") / 2.0
        perSegment = round(layers[0].getDimension("mult") / host.getDimension("mult") * segment)
        radii, layerMaterials, matrix = self._materials(host)
        key = (
            round(radius, 9),
            round(segment, 9),
            segments,
            perSegment,
            radii,
            tuple(id(m) for m in layerMaterials),
            id(matrix),
        )
        if key not in self._universes:
            particle = self._particleUniverse(radii, layerMaterials)
            centers = self._pack(radii[-1], radius, segment, perSegment)
            self._universes[key] = self._stack(
                particle, centers, radii[-1], radius, segment, segments, matrix
            )
            fraction = perSegment * _sphere(radii[-1]) / (math.pi * radius**2 * segment)
            runLog.info(
                f"OpenMC: packed {perSegment} particles per {segment:.3f} cm segment of "
                f"{host.name} (packing fraction {fraction:.4f}), stacked {segments} high"
            )
        return self._universes[key]

    def _fillSphere(self, host):
        """A universe of the particles packed in a spherical host, centered on the origin."""
        radius = host.getDimension("od") / 2.0
        count = round(particleLayers(host)[0].getParticleCount())
        radii, layerMaterials, matrix = self._materials(host)
        key = ("sphere", round(radius, 9), count, radii, tuple(map(id, layerMaterials)), id(matrix))
        if key not in self._universes:
            particle = self._particleUniverse(radii, layerMaterials)
            centers = self._pack(radii[-1], radius, None, count)
            pitch = self.options.particleMeshPitch or _autoPitch(radii[-1], count, _sphere(radius))
            n = max(1, math.ceil(2.0 * radius / pitch))
            trisos = [openmc.model.TRISO(radii[-1], particle, center=c) for c in centers]
            mesh = openmc.model.create_triso_lattice(
                trisos, (-radius,) * 3, (2.0 * radius / n,) * 3, (n, n, n), matrix
            )
            self._universes[key] = openmc.Universe(
                name="particles", cells=[openmc.Cell(name="packed sphere", fill=mesh)]
            )
            fraction = count * _sphere(radii[-1]) / _sphere(radius)
            runLog.info(
                f"OpenMC: packed {count} particles in {host.name} "
                f"(packing fraction {fraction:.4f}) on a {n}^3 mesh"
            )
        return self._universes[key]

    # ----------------------------------------------------------------------------------------
    def _particleUniverse(self, radii, layerMaterials):
        key = (radii, tuple(id(m) for m in layerMaterials))
        if key not in self._particleUniverses:
            spheres = [openmc.Sphere(r=r) for r in radii]
            cells = []
            for k, (sphere, mat) in enumerate(zip(spheres, layerMaterials, strict=True)):
                region = -sphere if k == 0 else (-sphere & +spheres[k - 1])
                cells.append(
                    openmc.Cell(name=mat.name if mat else "particle void", fill=mat, region=region)
                )
            # Outside the outermost sphere is the TRISO cell's container, never this universe.
            cells.append(openmc.Cell(name="outside particle", region=+spheres[-1]))
            self._particleUniverses[key] = openmc.Universe(name="particle", cells=cells)
        return self._particleUniverses[key]

    def _pack(self, particleRadius, radius, segment, count):
        """
        Particle centers in a cylinder of ``radius``, z from -segment/2 to segment/2, or in a
        sphere of ``radius`` if ``segment`` is None; cached.
        """
        seed = int(self.options.particleSeed)
        shape = "sphere" if segment is None else f"{segment:.9f}"
        tag = f"{particleRadius:.9f}-{radius:.9f}-{shape}-{count}-{seed}"
        cacheDir = os.path.join(self.options.runDir, "particle-packings")
        path = os.path.join(cacheDir, hashlib.sha1(tag.encode()).hexdigest()[:16] + ".npy")
        if os.path.exists(path):
            centers = np.load(path)
            if len(centers) == count:
                return centers
        if segment is None:
            region = -openmc.Sphere(r=radius)
        else:
            region = (
                -openmc.ZCylinder(r=radius)
                & +openmc.ZPlane(-segment / 2.0)
                & -openmc.ZPlane(segment / 2.0)
            )
        centers = openmc.model.pack_spheres(particleRadius, region, num_spheres=count, seed=seed)
        os.makedirs(cacheDir, exist_ok=True)
        np.save(path, centers)
        return centers

    def _stack(self, particle, centers, particleRadius, radius, segment, segments, matrix):
        """The packed segment, on its Cartesian mesh, stacked ``segments`` high from z = 0."""
        pitch = self.options.particleMeshPitch or _autoPitch(
            particleRadius, len(centers), math.pi * radius**2 * segment
        )
        nxy = max(1, math.ceil(2.0 * radius / pitch))
        nz = max(1, math.ceil(segment / pitch))
        pxy, pz = 2.0 * radius / nxy, segment / nz
        trisos = [openmc.model.TRISO(particleRadius, particle, center=c) for c in centers]
        mesh = openmc.model.create_triso_lattice(
            trisos, (-radius, -radius, -segment / 2.0), (pxy, pxy, pz), (nxy, nxy, nz), matrix
        )
        packed = openmc.Universe(
            name="packed segment", cells=[openmc.Cell(name="packed segment", fill=mesh)]
        )
        stack = openmc.RectLattice(name="segments")
        stack.lower_left = (-radius, -radius, 0.0)
        stack.pitch = (2.0 * radius, 2.0 * radius, segment)
        stack.universes = np.full((segments, 1, 1), packed, dtype=object)
        stack.outer = openmc.Universe(cells=[openmc.Cell(name="matrix", fill=matrix)])
        return openmc.Universe(
            name="particles", cells=[openmc.Cell(name="particle stack", fill=stack)]
        )


def _sphere(r):
    return 4.0 / 3.0 * math.pi * r**3


def _autoPitch(particleRadius, count, volume):
    """A mesh pitch that puts about ``PARTICLES_PER_ELEMENT`` particle centers in each element."""
    density = count / volume
    return max(2.0 * particleRadius, (PARTICLES_PER_ELEMENT / density) ** (1.0 / 3.0))
