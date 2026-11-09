Modeling the reactor
====================

Everything the plugin draws comes from the ARMI reactor model, at the hot dimensions ARMI uses
for its areas and volumes. Number densities therefore carry over unchanged: the volume of each
OpenMC cell is the volume ARMI holds the component's atoms in.

Blocks
------

A block becomes an OpenMC universe, infinite in z and in the plane; whatever places it clips it
to the block's footprint. ``openmcBlockModel`` picks the representation, and
``openmcHeterogeneousBlockFlags`` limits detail to some block types.

**Homogenized.** One cell filled with the block's volume-averaged number densities at its
volume-averaged temperature.

**Heterogeneous.** The components are laid out as ARMI describes them:

* **Pins** are components on the block's spatial grid (a ``MultiIndexLocation``). The components
  at each grid position, sorted by size, become the concentric layers of a pin universe (fuel,
  gap, clad; or moderator inside a guide tube). Pin universes fill an OpenMC lattice with the
  block grid's pitch and orientation: a ``RectLattice`` for Cartesian lattice maps, a
  ``HexLattice`` for hexagonal ones, including the pin grids ARMI creates for blocks given only
  a multiplicity such as ``mult: 217``.
* **Shells** are single centered shapes (circles, hexagons, squares, rectangles): ducts, inner
  ducts, cans, the intercoolant hexagon. They are nested by size around the pins.
* **Background.** Everything else fills the space that is left: ``DerivedShape`` coolant, which
  ARMI already defines that way, and the components OpenMC cannot draw exactly (helical wire
  wraps, for one), smeared in by area, which is what they physically displace.

Neighboring cells share their bounding surface object, so no particle is lost between two
coincident surfaces. Containment is checked flat by flat, which matters for a corners-up pin
bundle inside a flats-up duct, as in FFTF. A block that cannot be laid out consistently
(overlapping layers, pins wider than their pitch, pins inside a solid shell, a negative-area
component) is homogenized, with a warning that says why.

Materials
---------

Each component of each block gets its own OpenMC material, so a later change to one block
(burnup, temperature) never touches another. ARMI nuclide names map to the nuclides of the data
library: ``AM242M`` becomes ``Am242_m1``, and natural elements such as ``FE`` expand over the
isotopes the library has. Lumped fission products and dummy nuclides have no continuous-energy
data and are left out with a warning. Material temperatures are the component temperatures, and
``openmcTemperatureMethod`` decides how OpenMC treats temperatures between the library's.

S(α,β) tables are added from what a material is made of: water (H with O), heavy water,
graphite, beryllium, BeO and zirconium hydride. Homogenized LWR blocks therefore get water
thermal scattering too.

Particle fuel (TRISO)
---------------------

ARMI describes particle fuel with a ``particles:`` blueprint section and a host component that
holds the particles at a packing fraction (see the ARMI input documentation, *Particle fuel*).
Each particle layer is a ``ParticleLayer`` component of the block; the host is the matrix.
``openmcParticleModel`` decides how OpenMC sees them:

``explicit`` (the default)
   Every particle is drawn. Particle centers are randomly packed in the host with
   ``openmc.model.pack_spheres`` at ARMI's particle count, the particles are sorted onto a
   Cartesian mesh with ``openmc.model.create_triso_lattice`` (a neutron then checks only the few
   particles in its own mesh element), and the matrix fills the rest.

   A whole block height holds 10\ :sup:`5` and more particles per compact, too many to pack
   at once, so the host is packed in segments of about ``openmcParticleSegmentLength``
   (5 cm by default) that are stacked up the block, each a closed packing domain like the
   compacts stacked in a real fuel hole. Packings are cached under ``openmcRunDir`` by their
   parameters; ``openmcParticleSeed`` changes the packing and ``openmcParticleMeshPitch`` the
   mesh (by default about eight particle centers per mesh element).

   Hosts of identical dimensions, particle count and composition share one packed universe, so
   a core of fresh fuel blocks packs each kind of compact once: the MHTGR-350 core, 660 fuel
   blocks with 210 compacts each, packs one fuel and one burnable-poison compact and runs in
   under 500 MB.

``homogenized``
   Matrix and particles are smeared into one material by volume, compact by compact.

Explicit packing needs a solid circular host (a compact) or a solid sphere (a pebble's fuel
zone, packed whole); other hosts are smeared with a warning.

**Pebbles.** A ``SphericalBlock`` (a block whose largest component is a ``Sphere``) is a pebble's
unit cell: its spheres are drawn as concentric cells, and in XS mode the outer sphere gets a
white boundary, the usual closure of a pebble cell. Spherical cells have no place in a core
lattice and are refused there. Off-center components (``origin`` in the blueprints, such as a reserve shutdown
channel beside the pin lattice) are carved out of the pin region, and must not overlap a pin.

Cores
-----

The core's spatial grid becomes an OpenMC lattice of assembly universes; each assembly stacks
its blocks between the elevations ARMI holds, with void above and below. Assemblies keep their
own axial meshes, so cores with ``detailedAxialExpansion`` need no uniform-mesh conversion.
ARMI's flats-up hex grid (``cornersUp=False``) is OpenMC's ``'y'`` hex lattice; positions are
mapped through ``Lattice.find_element`` rather than by index arithmetic, so either orientation
and any Cartesian offset place correctly.

Symmetry:

* a full core is modeled as is;
* a quarter Cartesian core with reflective symmetry stays a quarter, cut by reflective planes
  through the origin (``openmcSymmetry: native``);
* any other partial core (third-core hex, periodic or eighth Cartesian) is unfolded to full core,
  placing each assembly's universe again at its symmetric positions.

A hexagonal core is bounded by a cylinder just around its assemblies, or of
``openmcCoreRadius``; a Cartesian core by the box of its lattice, or by a cylinder of
``openmcCoreRadius`` around it (assemblies in a pool or a drum). ``openmcOuterFill`` names a block
type whose composition fills the space between the assemblies and that boundary (a permanent side
reflector, the water of the pool). The block can be in the core or only in an assembly design of
the blueprints that the core map does not load.

Boundaries are vacuum by default; ``openmcTopBoundary``, ``openmcBottomBoundary`` and
``openmcRadialBoundary`` change them (C5G7 uses a reflective bottom). ``openmcAxialSlice``
builds a 2D model of the blocks at one elevation instead.

The initial fission source is uniform over the box around the fissionable blocks, not the
whole core, so shields and plena do not starve the first batch.

Results on the reactor
----------------------

These are mapped back after every run:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Parameter
     - Meaning
   * - ``core.keff``, ``core.openmcKeffStdDev``
     - Combined eigenvalue and its one-sigma uncertainty.
   * - block ``power``, ``pdens``
     - Fission energy deposited in the block (``kappa-fission``), W and W/cm³.
   * - block ``mgFlux``, ``flux``
     - Volume-integrated multigroup flux (n·cm/s, ARMI's fast-first group order) and the
       volume-averaged total flux (n/cm²/s).
   * - block ``linPowByPin``, ``openmcPinPowers``, ``openmcPinIndices``
     - With ``openmcPinPowers``: the power of every fissionable pin (W/cm, ordered like
       ``Block.getPinLocations``; and W with its block grid indices).

Tallies are per source neutron. They become absolute rates with the source rate that makes the
tallied fission energy equal the reactor power:

.. math::

   S = \frac{P}{e \sum_\text{model} \kappa\Sigma_f\phi}

Each block's tally is divided by the number of copies of the block in the model (an unfolded
core holds most assemblies several times) and, for an unfolded core, by the block's symmetry
factor, since ARMI stores only its share of a block cut by a symmetry line. In a reflective cut
the tallied block already is that share.
