# FFTF with OpenMC

The FFTF isothermal critical configuration (`FFTF.yaml`, `FFTF-blueprints.yaml`, from
terrapower/fftf-isothermal-model, Apache 2.0) run with the armi-openmc plugin. The benchmark is
FFTF-LMFR-RESR-001 (Bess et al., INL/EXT-09-16524, March 2010), whose benchmark-model k-eff is
0.9993 ± 0.0021 at 478 K (§3.1.5); its MCNP5/ENDF-B-VII.0 sample calculation gives
1.0031 ± 0.0001. The pin-detailed run here is (C−E)/E = −0.35 ± 0.21 %, with different data and
the ARMI inputs at 453 K.

| Settings | Blocks in detail | k (40k × 200, ENDF/B-VIII.1, 453 K) | Δρ |
|---|---|---|---|
| `FFTF-openmc-homog.yaml` | none | 0.98470 ± 0.00034 | |
| `FFTF-openmc-fuelonly.yaml` | fuel (217-pin bundles) | 0.99235 ± 0.00030 | +783 pcm from fuel heterogeneity |
| `FFTF-openmc.yaml` | fuel and absorbers (61 B4C pins inside an inner duct) | 0.99576 ± 0.00032 | +345 pcm more from absorber self-shielding |

Wire wraps are smeared into the sodium in every case. Blocks that the blueprints already give
as one homogeneous hexagon (shields, reflectors, plena, test assemblies) stay that way. The core
is full (313 assemblies; 2097 blocks; 6257 materials with every block detailed) and keeps each
assembly's own axial blocks.

The pin-detailed run took about 17 minutes on 12 threads, sharing the machine with other runs.

```bash
python -m armi_openmc run FFTF-openmc.yaml       # through ARMI; block power/flux in FFTF-openmc.h5
python -m armi_openmc openmc-core FFTF-openmc-homog.yaml --label fftf-homog
python -m armi_openmc openmc-core FFTF-openmc-fuelonly.yaml --label fftf-fuelonly
```

`python doc/make_figures.py` (from the repository root) draws the assembly power map from
`FFTF-openmc.h5`; see the Benchmarks page of the docs.
