Cross sections
==============

With ``xsKernel: OpenMC`` the plugin is ARMI's lattice physics code. For every representative
block of the cross section group manager it builds an OpenMC model, runs it, reduces the tallies
to microscopic group constants, and gathers everything into one ISOTXS.

Models
------

**Infinite lattice.** The representative block alone, its footprint with reflective sides.

**Colorset.** When the group's cross section settings name a ``driverID`` (and
``openmcXSColorset`` is on), the target sits among driver blocks:

* Cartesian: a 3×3 lattice with the target in the middle and reflective sides, which is an exact
  infinite checkerboard of the colorset;
* hexagonal: the target with one ring of drivers inside a cylinder of smeared driver with a white
  boundary. Hexagonal colorsets do not tile, so this approximates an infinite driver region.

A block with no fissionable material and no driver runs as a fixed-source problem with a
fission-spectrum source, which gives it a slowing-down spectrum.

Nuclides ARMI expects in the library but that the block lacks are added to every material at
``openmcXSTraceDensity``. They change nothing in the transport but are seen by the tallies.

Reduction
---------

For every nuclide :math:`i` the tallies give its reaction rates :math:`R_{i,g}` over the target
block, and the flux tally gives the volume-integrated group flux :math:`\Phi_g`. The
flux-volume weighted microscopic cross section is the one that reproduces the block's reaction
rates when multiplied by its block-averaged number density :math:`\bar N_i`:

.. math::

   \sigma_{i,g} = \frac{R_{i,g}}{\bar N_i \, \Phi_g}

* **Reaction cross sections** (total, fission, capture, (n,2n), (n,α), (n,p), (n,d), (n,t)),
  ν, the fission energy and the group velocities come from collision (expected-value)
  estimators. They score every nuclide at every collision, so trace nuclides converge with the
  flux, and they avoid the per-nuclide lookups a track-length tally would do on every short
  track segment of a pin lattice; for the FFTF fuel block that is about six times faster.
* **Scattering** takes its magnitude from the expected scattering rate and its shape (where
  neutrons go, and how many come out) from an analog energy-out tally. The matrix is
  production-weighted, (n,2n) included, and is written as the single P0 scattering block. ARMI
  counts (n,2n) as absorption, so the neutron balance closes. A nuclide too rare to have analog
  collisions keeps its scattering in-group.
* **Transport** cross sections subtract the analog P1 scattering moment from the total.
* **Fission spectra** come from analog outgoing-energy tallies; a nuclide without analog
  fissions takes the block's spectrum.
* **Empty groups.** A fast spectrum never reaches the thermal groups of a broad structure. Those
  groups take the constants of the nearest group the spectrum reaches, with in-group
  scattering, so no group of the library is singular.

Elements are tallied isotope by isotope and recombined by abundance. Each library also gets a
``<library>.flux.json`` beside it, with the group flux each XS ID was weighted with; this is
the weighting ``compare-xs`` collapses with.

Checks
------

Two checks run on every change (see :doc:`benchmarks` for the numbers):

* **Self-consistency.** With flux-volume weighted constants an infinite lattice keeps its neutron
  balance exactly. The k-infinity of the homogenized medium rebuilt from the written ISOTXS,
  through ARMI's own macroscopic cross section code, must therefore match OpenMC's k-infinity
  of the block within statistics. That tests densities, normalization, scattering orientation
  and band storage, fission spectra and energy order in one number (``armi_openmc.verify``,
  ``test_solver.py``).
* **An independent code.** ARMI's MC2-3 test library was made from a mixture whose MC2-3 input
  ships with ARMI; the plugin makes the same library and compares nuclide by nuclide.
