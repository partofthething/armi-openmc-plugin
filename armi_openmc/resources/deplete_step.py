"""
Deplete one ARMI time step with openmc.deplete.

Runs where OpenMC's shared library is (the OpenMC container, or a native install), in the run
directory that holds ``model.xml`` and ``depletion.json``. Written by
:py:mod:`armi_openmc.depletion` and copied there; it imports nothing from ARMI or the plugin.

``depletion.json``::

    {"chain": path, "integrator": "CECMIntegrator", "days": float, "power": W,
     "reactionRates": {"mode": "flux" | "direct", "groups": "SHEM-361",
                       "nuclides": [...], "reactions": [...]},
     "materials": {"<id>": {"volume": cm3, "nuclides": {"U235": atoms/b-cm, ...}}}}

In the flux mode, reaction rates come from collapsing each nuclide's cross sections with the group
flux, except for the listed nuclides and reactions, which are tallied directly.

The materials listed are the depletable ones. Their full inventory, including nuclides that have
no cross section data (OpenMC carries those through decay only), replaces what ``model.xml`` says.
Writes ``depletion-result.json``: keff at the start of the step and every depletable material's
number densities (atoms/b-cm) at its end.
"""

import json
import sys

import openmc
import openmc.deplete
import openmc.mgxs


def main(specPath="depletion.json"):
    with open(specPath) as f:
        spec = json.load(f)
    model = openmc.Model.from_model_xml("model.xml")
    burnable = spec["materials"]
    for mat in model.materials:
        entry = burnable.get(str(mat.id))
        if entry is None:
            # OpenMC makes any material with actinides depletable; only the listed ones are
            mat.depletable = False
            continue
        for nuc in list(mat.get_nuclides()):
            mat.remove_nuclide(nuc)
        for nuc, density in entry["nuclides"].items():
            mat.add_nuclide(nuc, density, "ao")
        mat.set_density("sum")
        mat.depletable = True
        mat.volume = entry["volume"]

    rates = spec.get("reactionRates", {"mode": "direct"})
    kwargs = {}
    if rates["mode"] == "flux":
        kwargs["reaction_rate_mode"] = "flux"
        kwargs["reaction_rate_opts"] = {
            "energies": openmc.mgxs.GROUP_STRUCTURES[rates["groups"]],
            "nuclides": rates["nuclides"],
            "reactions": rates["reactions"],
        }
    operator = openmc.deplete.CoupledOperator(
        model,
        spec["chain"],
        normalization_mode="fission-q",
        fission_yield_mode="constant",
        **kwargs,
    )
    integrator = getattr(openmc.deplete, spec["integrator"])(
        operator, [spec["days"]], power=[spec["power"]], timestep_units="d"
    )
    integrator.integrate(final_step=False, path="depletion_results.h5")

    results = openmc.deplete.Results("depletion_results.h5")
    _, keff = results.get_keff()
    end = results[-1]
    out = {"keff": [float(keff[0][0]), float(keff[0][1])], "materials": {}}
    for matId in burnable:
        volume = end.volume[matId]
        # atoms by (material, nuclide); older OpenMC also had a leading stage axis
        data = end.data if end.data.ndim == 2 else end.data[0]
        row = data[end.index_mat[matId]]
        densities = {}
        for nuc, j in end.index_nuc.items():
            atoms = float(row[j])
            if atoms > 0.0:
                densities[nuc] = atoms / volume * 1e-24
        out["materials"][matId] = densities
    with open("depletion-result.json", "w") as f:
        json.dump(out, f)


if __name__ == "__main__":
    main(*sys.argv[1:])
