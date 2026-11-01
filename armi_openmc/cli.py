"""
Command line entry points, available as ``python -m armi_openmc <command>``.

``openmc-core``
    Build the whole-core OpenMC model of a case, run it and report keff and block powers.
``openmc-xs``
    Make an ISOTXS library with OpenMC for the case's cross section groups, and optionally
    compare it with another ISOTXS.
``compare-xs``
    Compare two ISOTXS libraries nuclide by nuclide.
"""

import json
import os

import numpy as np
from armi import runLog
from armi.cli.entryPoint import EntryPoint


class BuildOpenMCModel(EntryPoint):
    """Build (and by default run) the whole-core OpenMC model of a case."""

    name = "openmc-core"
    settingsArgument = "required"

    def addOptions(self):
        self.parser.add_argument("--no-run", action="store_true", help="Only write the model.")
        self.parser.add_argument("--plot", action="store_true", help="Also write geometry plots.")
        self.parser.add_argument(
            "--label", default="core", help="Name of the run folder under openmcRunDir."
        )

    def invoke(self):
        from armi import init

        from armi_openmc import runner
        from armi_openmc.core import CoreModelBuilder
        from armi_openmc.options import OpenMCOptions
        from armi_openmc.results import CoreResults

        o = init(cs=self.cs)
        r = o.r
        opts = OpenMCOptions.fromSettings(
            self.cs, plotGeometry=self.args.plot or self.cs["openmcPlotGeometry"]
        )
        workdir = os.path.join(opts.runDir, self.args.label)
        os.makedirs(workdir, exist_ok=True)
        coreModel = CoreModelBuilder(r.core, opts).build()
        coreModel.model.export_to_model_xml(os.path.join(workdir, "model.xml"))
        if opts.plotGeometry:
            runner.run(opts, workdir, args=["-p"], logName="plot.out")
        if self.args.no_run:
            runLog.important(f"OpenMC model written to {workdir}")
            return 0
        results = CoreResults(coreModel, runner.run(opts, workdir))
        results.apply(r, self.cs["power"])
        summary = coreSummary(r, results)
        with open(os.path.join(workdir, "results.json"), "w") as f:
            json.dump(summary, f, indent=1)
        runLog.important(
            f"keff = {results.keff:.5f} +/- {results.keffStd:.5f}; "
            f"results in {workdir}/results.json"
        )
        return 0


def coreSummary(r, results):
    """keff and assembly/block powers as plain data."""
    assemblies = []
    for a in r.core:
        assemblies.append(
            {
                "name": a.getName(),
                "location": str(a.getLocation()),
                "type": a.getType(),
                "power": float(sum(b.p.power for b in a)),
                "blocks": [_blockSummary(b) for b in a],
            }
        )
    return {"keff": results.keff, "keffStdDev": results.keffStd, "assemblies": assemblies}


def _blockSummary(b):
    out = {
        "name": b.getName(),
        "type": b.getType(),
        "z": [float(b.p.zbottom), float(b.p.ztop)],
        "power": float(b.p.power),
        "flux": float(b.p.flux),
    }
    if b.p.openmcPinPowers is not None:
        out["pinPowers"] = [float(v) for v in b.p.openmcPinPowers]
        out["pinIndices"] = [[int(i) for i in ij] for ij in b.p.openmcPinIndices]
    return out


class GenerateXS(EntryPoint):
    """Make the case's ISOTXS library with OpenMC."""

    name = "openmc-xs"
    settingsArgument = "required"

    def addOptions(self):
        self.parser.add_argument("--output", default="ISOTXS-openmc", help="Library to write.")
        self.parser.add_argument("--compare", help="An ISOTXS to compare the new library with.")
        self.parser.add_argument("--xs-ids", nargs="+", help="Only these cross section IDs.")

    def invoke(self):
        from armi import init

        from armi_openmc import latticePhysics, mgxs

        o = init(cs=self.cs)
        xsGroups = o.getInterface("xsGroups")
        xsGroups.createRepresentativeBlocks()
        blocks = list(xsGroups.representativeBlocks.values())
        if self.args.xs_ids:
            blocks = [b for b in blocks if b.getMicroSuffix() in self.args.xs_ids]
        lib = latticePhysics.generateLibrary(o.r, self.cs, blocks)
        mgxs.writeLibrary(lib, self.args.output)
        if self.args.compare:
            compareLibraries(self.args.compare, self.args.output)
        return 0


class CompareXS(EntryPoint):
    """Compare two ISOTXS libraries nuclide by nuclide."""

    name = "compare-xs"

    def addOptions(self):
        self.parser.add_argument("reference", help="Reference ISOTXS.")
        self.parser.add_argument("other", help="ISOTXS to compare with the reference.")
        self.parser.add_argument(
            "--nuclides", nargs="+", help="Only these nuclide labels (e.g. U235AA)."
        )
        self.parser.add_argument(
            "--flux", help="Group flux JSON to collapse with (default: <other>.flux.json)."
        )

    def invoke(self):
        compareLibraries(self.args.reference, self.args.other, self.args.nuclides, self.args.flux)
        return 0


def compareLibraries(referencePath, otherPath, nuclides=None, fluxPath=None, out=print):
    """
    Print one-group and group-wise differences of the main cross sections of two ISOTXS files.

    Nuclides are matched by label (e.g. ``U235AA``). One-group values collapse with the flux the
    OpenMC library was made with (``<library>.flux.json``, written beside it), which is the
    weighting that means something; without one, flat in lethargy. The worst group difference
    looks only at groups holding at least 1% of the flux.
    """
    from armi.nuclearDataIO.cccc import isotxs

    ref = isotxs.readBinary(referencePath)
    other = isotxs.readBinary(otherPath)
    candidates = [fluxPath, f"{otherPath}.flux.json", f"{referencePath}.flux.json"]
    flux = None
    for path in candidates:
        if path and os.path.exists(path):
            with open(path) as f:
                flux = {k: np.asarray(v) for k, v in json.load(f).items()}
            out(f"Collapsing with the group flux in {path}")
            break
    rows = libraryDifferences(ref, other, nuclides, flux)
    out(
        f"{'nuclide':10s} {'reaction':12s} {'ref 1g':>11s} {'other 1g':>11s} "
        f"{'1g diff %':>9s} {'max |grp| %':>11s}"
    )
    for row in rows:
        out(
            f"{row['nuclide']:10s} {row['reaction']:12s} {row['ref']:11.4e} {row['other']:11.4e} "
            f"{row['diff']:9.2f} {row['maxGroupDiff']:11.2f}"
        )
    return rows


REACTIONS = {
    "total": lambda m: np.ravel(m.total[:, 0]) if np.ndim(m.total) == 2 else m.total,
    "transport": lambda m: (
        np.ravel(m.transport[:, 0]) if np.ndim(m.transport) == 2 else m.transport
    ),
    "fission": lambda m: m.fission,
    "nu": lambda m: m.neutronsPerFission,
    "capture": lambda m: m.nGamma,
    "n2n": lambda m: m.n2n,
    "scatter": lambda m: np.asarray(m.getTotalScatterMatrix().sum(axis=0)).ravel(),
}


def libraryDifferences(ref, other, nuclides=None, flux=None):
    """
    Rows of (nuclide, reaction, one-group values and differences) for the common nuclides.

    Parameters
    ----------
    flux : dict, optional
        XS ID -> group flux (fast group first) to collapse with. Flat in lethargy without it.
    """
    upper = np.asarray(ref.neutronEnergyUpperBounds)
    lower = np.append(upper[1:], 1e-5)
    lethargy = np.log(upper / lower)
    labels = [n for n in ref.nuclideLabels if n in other.nuclideLabels]
    if nuclides:
        labels = [n for n in labels if n in nuclides]
    rows = []
    for label in labels:
        xsID = label[-2:]
        weight = np.asarray(flux[xsID], float) if flux and xsID in flux else lethargy.copy()
        weight = weight / weight.sum()
        significant = weight >= 0.01
        a, b = ref[label].micros, other[label].micros
        for reaction, get in REACTIONS.items():
            try:
                va, vb = np.asarray(get(a), float), np.asarray(get(b), float)
            except (AttributeError, TypeError, ValueError):
                continue
            if va.shape != vb.shape or not np.any(va):
                continue
            if reaction == "nu":  # fission-weighted
                fa, fb = np.asarray(a.fission, float), np.asarray(b.fission, float)
                oneA = float(weight @ (va * fa) / (weight @ fa)) if weight @ fa else 0.0
                oneB = float(weight @ (vb * fb) / (weight @ fb)) if weight @ fb else 0.0
            else:
                oneA, oneB = float(weight @ va), float(weight @ vb)
            rel = np.divide(vb - va, va, out=np.zeros_like(va), where=np.abs(va) > 0)
            rows.append(
                {
                    "nuclide": label,
                    "reaction": reaction,
                    "ref": oneA,
                    "other": oneB,
                    "diff": 100.0 * (oneB - oneA) / oneA if oneA else 0.0,
                    "maxGroupDiff": 100.0 * float(np.max(np.abs(rel[significant]), initial=0.0)),
                }
            )
    return rows
