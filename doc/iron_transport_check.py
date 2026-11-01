"""
Is the Fe-56 transport difference between OpenMC and MC2-3 the flux-vs-current weighting?

MC2-3 writes a current-weighted (consistent-P1) transport cross section; the plugin's is
flux-weighted (total minus the P1 scattering moment). In an infinite medium with narrow
resonances the flux goes as 1/(E Sigma_t) and the current as 1/(E Sigma_t^2), so a nuclide with
deep windows in the mixture's total has a smaller current-weighted than flux-weighted total.
This computes that difference for Fe-56 from the pointwise ENDF/B-VIII.1 data of the MC2-3
fixture mixture, prints it next to the observed difference of the two libraries, and saves it
(``examples/mc2v3-fixture/iron-weighting.json``) for ``make_figures.py`` to draw.

    python doc/iron_transport_check.py
"""

import json
import os
import sys

import numpy as np
import openmc.data

sys.path.insert(0, os.path.dirname(__file__))
import make_figures  # noqa: E402

DATA = os.environ.get(
    "OPENMC_CROSS_SECTIONS", "/mnt/storage/nuclear_data/openmc/endfb-viii.1-hdf5/cross_sections.xml"
)
MIXTURE = {  # atoms/b-cm, from ARMI's mc2v3-AA.inp
    "U235": 1e-3, "U238": 1e-3, "Pu239": 1e-4, "Fe54": 1e-3, "Fe56": 1e-2, "Fe57": 1e-4,
    "Fe58": 1e-5, "Na23": 1e-3, "Zr90": 1e-3, "Zr91": 1e-4, "Zr92": 1e-4, "Zr94": 1e-4,
    "Zr96": 1e-5,
}  # fmt: skip


def totals(energy):
    library = openmc.data.DataLibrary.from_xml(DATA)
    out = {}
    for name in MIXTURE:
        nuc = openmc.data.IncidentNeutron.from_hdf5(library.get_by_material(name)["path"])
        temp = min(nuc.temperatures, key=lambda t: abs(int(t[:-1]) - 873))
        xs = nuc[1].xs[temp]
        out[name] = np.interp(energy, xs.x, xs.y)
    return out


def main():
    ref, omc, weight, upper, lower = make_figures._mc2v3Libraries()
    keep = np.flatnonzero(weight >= 1e-5)
    energy = np.geomspace(lower[keep].min(), upper[keep].max(), 400_000)
    sigma = totals(energy)
    macro = sum(MIXTURE[n] * sigma[n] for n in MIXTURE)
    flux = 1.0 / (energy * macro)  # narrow-resonance flux
    current = flux / macro  # and current
    trA = np.ravel(ref["FE56AA"].micros.transport[:, 0])
    trB = np.ravel(omc["FE56AA"].micros.transport[:, 0])
    print(
        f"{'group (keV)':>18s} {'flux-wtd tot':>12s} {'curr-wtd tot':>12s} "
        f"{'predicted':>9s} {'observed':>9s}  (barns, OpenMC - MC2-3 transport)"
    )
    predicted = {}
    for g in keep:
        inGroup = (energy >= lower[g]) & (energy < upper[g])
        e, s = energy[inGroup], sigma["Fe56"][inGroup]
        byFlux = np.trapezoid(s * flux[inGroup], e) / np.trapezoid(flux[inGroup], e)
        byCurrent = np.trapezoid(s * current[inGroup], e) / np.trapezoid(current[inGroup], e)
        predicted[int(g)] = float(byFlux - byCurrent)
        print(
            f"{lower[g] / 1e3:8.3g}-{upper[g] / 1e3:8.3g} {byFlux:12.3f} {byCurrent:12.3f} "
            f"{byFlux - byCurrent:9.3f} {trB[g] - trA[g]:9.3f}"
        )
    out = make_figures.ROOT / "examples" / "mc2v3-fixture" / "iron-weighting.json"
    out.write_text(json.dumps({"predictedTransportGapBarns": predicted}, indent=1))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
