"""
Compare a 2D OpenMC run of C5G7 with the benchmark's 2D reference.

Run the model first::

    python -m armi_openmc openmc-core c5g7-2d.yaml --label c5g7-2d

The reference (NEA/NSC/DOC(2003)16, MCNP with the benchmark's 7-group cross sections) is for the
multigroup problem. This model is continuous-energy ENDF/B-VIII.1 built from the C5
compositions the 7-group data were made from, so it checks the plugin's geometry and mapping, not
the benchmark's numbers to the pcm: expect agreement to a few percent, not to statistics.
"""

import json
import sys

import numpy as np

REFERENCE = {
    "keff": 1.18655,
    # assembly power, normalized so the average fuel pin is 1 (264 pins per assembly)
    "inner UO2": 492.8,
    "MOX": 211.7,
    "outer UO2": 139.1,
    "max pin": 2.498,
    "min pin": 0.232,
}


def main(path):
    data = json.load(open(path))
    fuel = []
    for a in data["assemblies"]:
        pins = [p for b in a["blocks"] for p in b.get("pinPowers", [])]
        if pins:
            fuel.append((a["location"], a["type"], np.array(pins)))
    allPins = np.concatenate([p for _, _, p in fuel])
    # The fission chamber's trace U-235 is a fissionable pin with ~no power; drop it.
    norm = allPins[allPins > 1e-3 * allPins.max()].mean()
    print(
        f"keff  OpenMC {data['keff']:.5f} +/- {data['keffStdDev']:.5f}   "
        f"reference {REFERENCE['keff']}"
    )
    for loc, typ, pins in sorted(fuel):
        pins = pins[pins > 1e-3 * pins.max()]
        print(
            f"{typ:6s} {loc:12s} assembly power {pins.sum() / norm:7.1f}  "
            f"max pin {pins.max() / norm:.3f}  min {pins.min() / norm:.3f}"
        )
    print("reference: inner UO2 492.8, MOX 211.7, outer UO2 139.1, max pin 2.498, min pin 0.232")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "openmc-runs/c5g7-2d/results.json")
