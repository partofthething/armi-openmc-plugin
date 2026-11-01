# armi-openmc

An [ARMI](https://partofthething.github.io/armi/) plugin that builds [OpenMC](https://docs.openmc.org)
Monte Carlo models straight from an ARMI reactor, runs them, and puts the results back on the
reactor. It is written against the partofthething fork of ARMI (2.0+).

Documentation: https://partofthething.github.io/armi-openmc-plugin/
Source and issues: https://github.com/partofthething/armi-openmc-plugin

It does two jobs:

| Mode | Setting | What it does |
|---|---|---|
| **Global flux** | `neutronicsKernel: OpenMC` | A whole-core eigenvalue run at every time node: `keff` (with its standard deviation), block `power`, `pdens`, `flux` and multigroup `mgFlux`, and optionally pin powers (`linPowByPin`, `openmcPinPowers`). |
| **Cross sections** | `xsKernel: OpenMC` | For each cross section group, a model of its representative block (an infinite lattice, or a colorset with its `driverID`), tallied into microscopic multigroup cross sections in the case's group structure and written as an **ISOTXS** that any ARMI flux solver can use. |

## Geometry

* **Cores:** hexagonal (flats-up or corners-up) and Cartesian grids, as an OpenMC lattice of
  assembly universes. Every assembly keeps its own block boundaries, so cores with
  `detailedAxialExpansion` or non-uniform meshes run as they are, with no axial remeshing.
* **Symmetry:** a quarter Cartesian core with reflective symmetry is modeled as the quarter, cut by
  reflective planes. Third-core hex and the other partial cores are unfolded to full core, and
  tallies are divided back onto ARMI's blocks, including the split center assembly.
* **Blocks**, chosen with `openmcBlockModel`:
  * `heterogeneous`: pins on the block's spatial grid (Cartesian lattice maps, or the hexagonal
    pin grids ARMI makes for `mult: 217`-style blocks) as concentric layers in a pin lattice.
    Ducts, cans and intercoolant are nested shells around the pins. `DerivedShape` coolant fills
    the rest, and wire wraps and other undrawable shapes are smeared into it by area.
    `openmcHeterogeneousBlockFlags` limits the detail to, say, `[fuel, control]`.
  * `homogenized`: one material per block.

  A block that cannot be laid out (overlapping layers, pins wider than the pitch) is homogenized
  with a warning saying why. Volumes match the ARMI component areas, so number densities carry
  over unchanged.
* **2D:** `openmcAxialSlice: <z>` models the blocks at one elevation, with reflective top and
  bottom.
* **Materials:** ARMI nuclides map to the library's (natural elements expand over the isotopes it
  has). Temperatures come from the components. S(α,β) tables are added automatically for water,
  heavy water, graphite, beryllium, BeO and ZrH.

## Running OpenMC

`openmcRunner: auto` uses an `openmc` executable on the `PATH`, or else the `openmc/openmc`
container: docker mounts the run directory and the nuclear data at their own paths and runs as you.
Point `openmcCrossSections` (or `OPENMC_CROSS_SECTIONS`) at a `cross_sections.xml`. Each run gets a
folder under `openmcRunDir` with `model.xml`, the statepoint, `openmc.out`, and geometry plots when
`openmcPlotGeometry` is on. OpenMC's batch output is echoed into the ARMI log as `[openmc] ...`
lines.

## Install

PyPI has no OpenMC wheels, so the Python API comes from its release tag; the solver can stay in
docker.

```bash
git clone https://github.com/partofthething/armi-openmc-plugin.git && cd armi-openmc-plugin
uv venv && uv pip install -e .          # armi from the fork's git main, openmc 0.15.3 from git
uv pip install --no-deps -e ../armi     # optional: a local, editable ARMI fork
```

Installing this package pulls ARMI from git, which replaces an editable ARMI that was already in
the environment. Reinstall the local one afterwards, as the second line does.

## Use

```bash
python -m armi_openmc run case.yaml            # an ARMI run with the plugin loaded
python -m armi_openmc gui case.yaml            # the input editor: an "OpenMC" settings section,
                                               # and Run (Ctrl+Enter) to run the case and follow its log
python -m armi_openmc openmc-core case.yaml    # just build/run the core model and write results.json
python -m armi_openmc openmc-xs case.yaml --compare ISOTXS-other   # make an ISOTXS and diff it
python -m armi_openmc compare-xs ISOTXS-a ISOTXS-b
```

Another app can register `armi_openmc.plugin.OpenMCPlugin` itself instead.

## Documentation

Published at https://partofthething.github.io/armi-openmc-plugin/ (built and deployed by
`.github/workflows/docs.yaml` on every push to main). Sphinx docs (Book theme) in `doc/`: usage, how the models and cross sections are made, the
benchmarks with figures, the settings reference, and the roadmap (depletion, variance reduction,
richer coupling back onto ARMI).

```bash
uv pip install -e ".[docs]"
python doc/make_figures.py                 # after running the examples
sphinx-build -b html doc doc/_build/html
```

## Examples and checks

| Case | Result |
|---|---|
| `examples/c5g7/c5g7-2d.yaml`: C5G7 2D, quarter core, pin by pin, 50k × 250 | k = 1.17095 ± 0.00032. Assembly powers 479 / 218 / 141, against the benchmark's MCNP multigroup reference 493 / 212 / 139 (inner UO2 / MOX / outer UO2); max pin 2.436 vs 2.498. The reference uses the 7-group benchmark data, while this run is continuous-energy ENDF/B-VIII.1 on the C5 compositions, so expect a few percent, not statistical agreement (`examples/c5g7/compare_2d.py`). |
| `examples/c5g7/c5g7-openmc.yaml`: C5G7 3D unrodded through `armi run` | k = 1.16637 ± 0.00072 (20k × 150). |
| `examples/fftf/FFTF-openmc.yaml`: FFTF isothermal critical, full hex core, all pin blocks detailed (2097 blocks, 6257 materials) | k = 0.99576 ± 0.00032 (40k × 200) against the benchmark's 0.9993 ± 0.0021: (C−E)/E = −0.35 ± 0.21 %. Fuel heterogeneity is worth +783 pcm and absorber self-shielding +345 pcm (homogenized: 0.98470). See `examples/fftf/README-openmc.md`. |
| `examples/mc2v3-fixture`: ARMI's MC2-3 ISOTXS fixture mixture, OpenMC XS mode | U-235/U-238/Pu-239/Na-23 fission, ν, capture and total within 0.1–2% of MC2-3 (ENDF/B-VII.0 vs VIII.1); `comparison.txt` |
| XS self-consistency: the k-infinity of the homogenized medium built from the written ISOTXS, against OpenMC's k-infinity of the block | C5G7 UO2 assembly (homogenized, 33 groups with upscatter): 1.25812 vs 1.2586 ± 0.0006, run as a test (`test_solver.py`). FFTF inner fuel block (217 pins in detail, fast spectrum, 5 empty thermal groups filled): 1.51092 vs 1.5121 ± 0.0009, in about 4 minutes on 6 threads. |

Tests: `pytest armi_openmc/tests`. Tests marked `openmc` run the solver at low statistics.

## How the cross sections are reduced

See the docstring of `armi_openmc/mgxs.py`. In short:
* Reaction cross sections are flux-volume weighted collision (expected-value) tallies divided by
  the block-averaged number density, so trace nuclides converge as well as the flux does.
* The scattering matrix (production-weighted, P0, stored as the elastic block) takes its
  magnitude from the expected scattering rate and its shape from an analog estimator. It is in-group for nuclides
  too rare to have analog statistics.
* Transport uses the analog P1 moment.
* χ is an analog outgoing-energy tally, falling back to the block's χ.

Nuclides missing from a block are added at `openmcXSTraceDensity` so that the library covers
everything ARMI expects. Groups the block's spectrum never reaches (the thermal groups of a fast
block) get the constants of the nearest group that it does reach, with in-group scattering, so
the library has no singular groups.

## Limitations

* No depletion through `openmc.deplete`: burnup stays ARMI's job, as with other kernels.
* Hexagonal colorsets are approximate: a ring of drivers inside a white cylinder of smeared driver.
* Only P0 scattering is written to ISOTXS (with transport-corrected totals); gamma libraries are not.
* Components with shapes other than circles, hexagons, squares and rectangles are smeared, not drawn.
