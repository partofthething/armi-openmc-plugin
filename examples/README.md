# Examples

All of these use ENDF/B-VIII.1 (`openmcCrossSections` points at it; change the path for your
machine) and OpenMC 0.15.3 in docker unless an `openmc` executable is on the `PATH`.

## C5G7 (`c5g7/`)

The C5G7 MOX benchmark as described in the ARMI fork's test inputs
(`armi/testing/reactors/c5g7`): a quarter core of UO2 and MOX 17×17 assemblies with a water
reflector, using the C5 compositions the benchmark's 7-group data were made from. Every assembly
is modeled pin by pin, and the quarter stays a quarter, cut by reflective planes.

```bash
cd c5g7
python -m armi_openmc openmc-core c5g7-2d.yaml --label c5g7-2d   # 2D: a slice through the fuel
python compare_2d.py                                              # against the benchmark reference
python -m armi_openmc run c5g7-openmc.yaml                        # 3D, through the ARMI interface stack
```

| Run | k | Assembly powers (inner UO2 / MOX / outer UO2) | Max pin |
|---|---|---|---|
| OpenMC 2D, 50k × 250 | 1.17095 ± 0.00032 | 479 / 218 / 141 | 2.436 |
| Benchmark reference (MCNP, 7-group) | 1.18655 | 492.8 / 211.7 / 139.1 | 2.498 |
| OpenMC 3D unrodded, reflective bottom, 20k × 150 | 1.16637 ± 0.00072 | | |

The reference solves the 7-group problem, while this is continuous energy, so a difference of a
few percent is the physics of the data, not the plugin. The two MOX assemblies, which sit at
symmetric positions, come out identical, as they should.

## FFTF (`fftf/`)

The FFTF isothermal critical configuration (control rods at their critical positions), copied
from [terrapower/fftf-isothermal-model](https://github.com/terrapower/fftf-isothermal-model)
(Apache 2.0, see `fftf/LICENSE.md`). The benchmark is FFTF-LMFR-RESR-001 (Bess et al.,
INL/EXT-09-16524), and its evaluated k-eff is in the IRPhEP handbook. The model is a full hex
core of 313 assemblies; fuel and absorber blocks are pin by pin (217-pin fuel bundles, 61-pin B4C
bundles inside an inner duct, wire wraps smeared into the sodium), and the rest are the
homogenized blocks the blueprints define.

```bash
cd fftf
python -m armi_openmc run FFTF-openmc.yaml                          # pin-detailed, writes FFTF-openmc.h5
python -m armi_openmc openmc-core FFTF-openmc-homog.yaml --label fftf-homog   # every block homogenized
python -m armi_openmc openmc-core FFTF-openmc-fuelonly.yaml --label fftf-fuelonly
```

| Run (40k × 200) | k |
|---|---|
| Every block homogenized | 0.98470 ± 0.00034 |
| Fuel pin-detailed | 0.99235 ± 0.00030 |
| Fuel and absorbers pin-detailed (2097 blocks, 6257 materials) | 0.99576 ± 0.00032 |

Details are in `fftf/README-openmc.md`.

The docs (`doc/benchmarks.rst`) compare these with the benchmark; `python doc/make_figures.py`
draws the assembly power map from `FFTF-openmc.h5`.

## MC2-3 fixture (`mc2v3-fixture/`)

ARMI's MC2-3 test library `mc2v3-AA.isotxs` was made from an infinite homogeneous fast mixture
whose MC2-3 input ships with ARMI. `compare_mc2v3.py` builds the same mixture as an ARMI block,
makes its ISOTXS with OpenMC (XS mode), and compares the two, collapsing with the OpenMC flux.
`comparison.txt` has the full table. For the main nuclides:

| | fission | ν | capture | total |
|---|---|---|---|---|
| U-235 | +0.3% | −0.1% | +1.8% | +1.5% |
| U-238 | +1.1% | −0.1% | −0.6% | 0.0% |
| Pu-239 | +0.2% | +0.2% | −1.7% | −1.3% |
| Na-23 | | | +0.3% | +0.3% |

These are ENDF/B-VIII.1 against MC2-3's ENDF/B-VII.0, with MC2-3 adding leakage through its
buckling search. Iron's transport cross section is 13–20% higher here, because MC2-3 weights
transport by current and this plugin by flux, which differ most in iron's resonance windows.
