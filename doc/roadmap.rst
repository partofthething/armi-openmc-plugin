Roadmap
=======

What the plugin does today is listed in :doc:`usage` and :doc:`modeling`. The plans, roughly in
order:

Depletion
---------

Activate depletion through ``openmc.deplete``, so the plugin can burn the reactor as well as
solve it:

* coupled transport-depletion over ARMI's cycles and burn steps, with OpenMC's integrators
  (predictor, CE/CM, CF4) as options;
* a depletion chain matched to the data library, chosen by setting;
* depleted number densities written back onto the ARMI components, with burnup, so ARMI's fuel
  management, database and reports see them, and the results restartable from the ARMI
  database;
* depletion at pin resolution where the model is pin by pin, and per block where it is
  homogenized.

Variance reduction
------------------

Activate OpenMC's variance reduction options, so large and deep-penetration problems converge in
reasonable time:

* weight windows, generated automatically (MAGIC, or FW-CADIS from a random-ray adjoint) for
  shielding, ex-core detector and vessel fluence problems;
* survival biasing and uniform fission site sampling, to flatten statistics across a core for
  pin powers;
* source convergence diagnostics and automatic inactive batch selection from the Shannon entropy;
* uncertainty targets: run until a keff or tally precision is reached instead of a fixed number of
  batches.

Bringing more results back onto ARMI for physics coupling
---------------------------------------------------------

Today the plugin maps keff, block power and power density, block-level multigroup flux and total
flux, and per-pin power onto the ARMI model. The aim is to make the ARMI model carry everything
another physics code needs from Monte Carlo:

* **power distributions** at every level, with their statistical uncertainties: pin-by-pin and
  axially resolved within blocks, plus heating from photons (coupled neutron-photon transport),
  for thermal-hydraulics and fuel performance coupling;
* **multigroup fluxes at pin level** and per component, for fluence, dpa and reaction rate edits;
* **number densities** from depletion, back on the components (above);
* temperature and density **feedback** from ARMI's thermal-hydraulics, through per-cell
  temperatures and densities that change without rebuilding the model, for tightly coupled
  iterations.

Also planned
------------

* MPI execution on clusters, next to the threaded and container runners.
* Current-weighted (consistent-P1) transport cross sections in XS mode, as MC2-3 writes them.
  The flux-weighted ones the plugin writes now overstate iron's transport cross section by up to
  tens of percent between 20 keV and 1 MeV (see :doc:`benchmarks`), which matters for diffusion
  coefficients in fast cores.
* Higher-order (P1+) scattering and gamma libraries (GAMISO/PMATRX) in XS mode, and better
  hexagonal colorsets.
* Shapes beyond circles, hexagons and rectangles drawn exactly instead of smeared.
* Pebble beds: pebbles packed in a vessel, with the same packing machinery as compacts. (A
  pebble's unit cell, with explicit particles, is done: see LEUPRO-1 in :doc:`benchmarks`.)
* Reactivity coefficient and control rod worth calculations driven from ARMI.
* A standing set of verification problems rerun on every change, in CI.
