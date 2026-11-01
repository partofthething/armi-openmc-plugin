"""
Make the benchmark figures for the docs from the example runs.

Run from the repository root after running the examples (see ``examples/README.md``)::

    python doc/make_figures.py

Writes ``doc/_static/c5g7-pin-power.png`` from the 2D C5G7 ``results.json``,
``doc/_static/fftf-assembly-power.png`` from the FFTF ARMI database, and the comparison with
ARMI's MC2-3 fixture library (``xs-mc2v3-*.png``) from ``examples/mc2v3-fixture``.
"""

import json
import math
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import RegularPolygon  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
STATIC = ROOT / "doc" / "_static"

# One hue, light to dark (sequential): magnitude only.
BLUES = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
SEQUENTIAL = LinearSegmentedColormap.from_list("seq-blue", BLUES)
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
EMPTY = "#f0efec"
GRID = "#e1e0d9"
SERIES = {"MC2-3 (ENDF/B-VII.0)": "#2a78d6", "OpenMC (ENDF/B-VIII.1)": "#eb6834"}
# Diverging: blue (OpenMC lower) - gray - red (OpenMC higher)
DIVERGING = LinearSegmentedColormap.from_list(
    "div", ["#184f95", "#6da7ec", "#f0efec", "#ee8a89", "#a62b2b"]
)


def _style(ax, fig):
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=8)


def _colorbar(fig, mappable, ax, label):
    bar = fig.colorbar(mappable, ax=ax, fraction=0.04, pad=0.02)
    bar.outline.set_visible(False)
    bar.ax.tick_params(colors=INK_2, labelsize=8)
    bar.set_label(label, color=INK_2, fontsize=9)
    return bar


def c5g7PinPower(results=ROOT / "examples/c5g7/openmc-runs/c5g7-2d/results.json"):
    data = json.loads(pathlib.Path(results).read_text())
    n = 17
    grid = np.full((3 * n, 3 * n), np.nan)
    for a in data["assemblies"]:
        ai, aj = (int(v) for v in a["location"].split("-")[:2])
        for b in a["blocks"]:
            for (i, j), p in zip(b.get("pinIndices", []), b.get("pinPowers", []), strict=True):
                grid[aj * n + j + 8, ai * n + i + 8] = p
    fuel = grid[np.isfinite(grid)]
    trace = 1e-3 * np.nanmax(grid)
    grid[grid < trace] = np.nan  # the fission chamber's trace U-235
    grid /= np.nanmean(grid[grid >= trace])

    fig, ax = plt.subplots(figsize=(6.4, 5.6), dpi=160)
    _style(ax, fig)
    ax.add_patch(plt.Rectangle((0, 0), 3 * n, 3 * n, color=EMPTY, zorder=0))
    image = ax.imshow(
        grid, origin="lower", cmap=SEQUENTIAL, extent=(0, 3 * n, 0, 3 * n), interpolation="nearest"
    )
    for k in range(1, 3):  # assembly boundaries
        ax.axhline(k * n, color=SURFACE, lw=2)
        ax.axvline(k * n, color=SURFACE, lw=2)
    peak = np.unravel_index(np.nanargmax(grid), grid.shape)
    ax.add_patch(
        plt.Rectangle((peak[1], peak[0]), 1, 1, fill=False, edgecolor=SURFACE, lw=1.6, zorder=3)
    )
    ax.text(
        0.5 * n,
        2.5 * n,
        f"peak pin {np.nanmax(grid):.3f}\n(outlined, inner UO2)",
        ha="center",
        va="center",
        color=INK_2,
        fontsize=8,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(
        f"C5G7 2D pin power, quarter core (k = {data['keff']:.5f})",
        color=INK,
        fontsize=10,
        loc="left",
    )
    for (cx, cy), name in {(0, 0): "UO2", (1, 0): "MOX", (0, 1): "MOX", (1, 1): "UO2"}.items():
        ax.text(
            (cx + 0.5) * n,
            (cy + 1) * n - 1.6,
            name,
            ha="center",
            color=INK,
            fontsize=8,
            bbox={"boxstyle": "round,pad=0.2", "fc": SURFACE, "ec": "none", "alpha": 0.85},
        )
    ax.text(2.5 * n, 1.5 * n, "water\nreflector", ha="center", va="center", color=MUTED, fontsize=8)
    _colorbar(fig, image, ax, "pin power / average fuel pin")
    out = STATIC / "c5g7-pin-power.png"
    fig.savefig(out, bbox_inches="tight", facecolor=SURFACE)
    print(f"wrote {out} ({fuel.size} pins)")


def fftfAssemblyPower(dbPath=ROOT / "examples/fftf/FFTF-openmc.h5"):
    import armi

    if not armi.isConfigured():
        from armi_openmc.app import OpenMCApp

        armi.configure(OpenMCApp())
    from armi.bookkeeping.db import databaseFactory

    with databaseFactory(str(dbPath), "r") as db:
        r = db.load(0, 0)
    pitch = r.core.spatialGrid.pitch
    radius = pitch / math.sqrt(3.0)
    centers, powers = [], []
    for a in r.core:
        x, y, _ = a.spatialLocator.getGlobalCoordinates()
        centers.append((x, y))
        powers.append(sum(b.p.power for b in a) / 1e6)
    powers = np.array(powers)
    top = powers.max()

    fig, ax = plt.subplots(figsize=(6.4, 6.0), dpi=160)
    _style(ax, fig)
    floor = 3.0  # just under the lowest fuel assembly, so the core's own variation shows
    norm = plt.Normalize(floor, top)
    for (x, y), p in zip(centers, powers, strict=True):
        # flat-topped assemblies: a vertex on the x axis
        face = SEQUENTIAL(norm(p)) if p > 0.01 else EMPTY
        ax.add_patch(
            RegularPolygon(
                (x, y),
                6,
                radius=radius,
                orientation=math.pi / 6,
                facecolor=face,
                edgecolor=SURFACE,
                lw=1.2,
            )
        )
        if p > 0.01:
            ax.text(
                x,
                y,
                f"{p:.1f}",
                ha="center",
                va="center",
                fontsize=4.6,
                color="#ffffff" if norm(p) > 0.45 else INK,
            )
    xs, ys = zip(*centers, strict=True)
    ax.set_xlim(min(xs) - pitch, max(xs) + pitch)
    ax.set_ylim(min(ys) - pitch, max(ys) + pitch)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(
        f"FFTF assembly power, MW (k = {r.core.p.keff:.5f})", color=INK, fontsize=10, loc="left"
    )
    mappable = plt.cm.ScalarMappable(norm=norm, cmap=SEQUENTIAL)
    _colorbar(fig, mappable, ax, "assembly power (MW); gray: no fuel")
    ax.text(
        0.0,
        -0.02,
        "Scale starts at 3 MW; the lowest fuel assembly makes 3.6 MW.",
        transform=ax.transAxes,
        color=MUTED,
        fontsize=7,
    )
    out = STATIC / "fftf-assembly-power.png"
    fig.savefig(out, bbox_inches="tight", facecolor=SURFACE)
    print(f"wrote {out}")


def _mc2v3Libraries():
    import armi

    if not armi.isConfigured():
        from armi_openmc.app import OpenMCApp

        armi.configure(OpenMCApp())
    from armi.nuclearDataIO.cccc import isotxs

    here = ROOT / "examples" / "mc2v3-fixture"
    fixture = pathlib.Path(armi.__file__).parent / "nuclearDataIO/tests/fixtures/mc2v3-AA.isotxs"
    ref = isotxs.readBinary(str(fixture))
    omc = isotxs.readBinary(str(here / "ISOAA-openmc"))
    flux = np.array(json.loads((here / "ISOAA-openmc.flux.json").read_text())["AA"])
    upper = np.asarray(ref.neutronEnergyUpperBounds)
    lower = np.append(upper[1:], 1e-5)
    return ref, omc, flux / flux.sum(), upper, lower


def _vector(micros, reaction):
    values = {
        "fission": micros.fission,
        "capture": micros.nGamma,
        "total": np.ravel(micros.total[:, 0]) if np.ndim(micros.total) == 2 else micros.total,
        "nu": micros.neutronsPerFission,
    }[reaction]
    return np.asarray(values, float)


def _minorLabel(value, _pos):
    """Label the 2 and 5 minor ticks of a log axis, for panels that span less than a decade."""
    if value <= 0:
        return ""
    mantissa = round(value / 10 ** math.floor(math.log10(value)), 6)
    return f"{value:g}" if mantissa in (2.0, 5.0) else ""


PANELS = [
    ("U235AA", "fission", "U-235 fission"),
    ("PU39AA", "fission", "Pu-239 fission"),
    ("U238AA", "capture", "U-238 capture"),
    ("U238AA", "total", "U-238 total"),
    ("FE56AA", "total", "Fe-56 total"),
    ("NA23AA", "total", "Na-23 total"),
]


def mc2v3CrossSections():
    """Group-wise cross sections of both libraries, and their difference, in small multiples."""
    ref, omc, weight, upper, lower = _mc2v3Libraries()
    # Groups that hold the spectrum (99.99% of the flux): below that the mixture has no neutrons
    # and both libraries are extrapolating.
    keep = weight >= 1e-5
    edges = np.append(upper[keep], lower[keep][-1])

    fig, axes = plt.subplots(
        4, 3, figsize=(9.6, 8.6), dpi=160, gridspec_kw={"height_ratios": [3, 1.3, 3, 1.3]}
    )
    fig.patch.set_facecolor(SURFACE)
    for k, (label, reaction, title) in enumerate(PANELS):
        row, col = 2 * (k // 3), k % 3
        top, bottom = axes[row, col], axes[row + 1, col]
        a, b = (
            _vector(ref[label].micros, reaction)[keep],
            _vector(omc[label].micros, reaction)[keep],
        )
        for (name, color), values in zip(SERIES.items(), (a, b), strict=True):
            top.stairs(values, edges, color=color, lw=2, label=name, baseline=None)
        diff = 100.0 * (b - a) / np.where(a > 0, a, np.nan)
        bottom.stairs(diff, edges, color=INK_2, lw=1.5, baseline=None)
        bottom.axhline(0.0, color="#c3c2b7", lw=1)
        oneA, oneB = float(weight[keep] @ a), float(weight[keep] @ b)
        top.set_title(
            f"{title}   1-group {100 * (oneB - oneA) / oneA:+.1f}%",
            color=INK,
            fontsize=9,
            loc="left",
        )
        for ax in (top, bottom):
            _style(ax, fig)
            ax.set_xscale("log")
            ax.set_xlim(edges[-1], edges[0])
            ax.grid(True, color=GRID, lw=0.6)
            ax.set_axisbelow(True)
        top.set_yscale("log")
        top.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        lo, hi = top.get_ylim()
        narrow = math.log10(hi / lo) < 1.3
        top.yaxis.set_minor_formatter(FuncFormatter(_minorLabel if narrow else lambda v, _: ""))
        top.tick_params(which="both", colors=MUTED, labelsize=8)
        top.tick_params(labelbottom=False)
        lim = max(5.0, np.nanpercentile(np.abs(diff), 95) * 1.2)
        bottom.set_ylim(-lim, lim)
        if col == 0:
            top.set_ylabel("barns", color=INK_2, fontsize=8)
            bottom.set_ylabel("diff %", color=INK_2, fontsize=8)
        if row == 2:
            bottom.set_xlabel("neutron energy (eV)", color=INK_2, fontsize=8)
    handles, names = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        names,
        loc="upper left",
        ncol=2,
        frameon=False,
        fontsize=9,
        bbox_to_anchor=(0.01, 1.0),
        labelcolor=INK,
    )
    fig.text(
        0.99,
        0.995,
        "diff % = (OpenMC - MC2-3) / MC2-3",
        ha="right",
        va="top",
        color=MUTED,
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = STATIC / "xs-mc2v3-groupwise.png"
    fig.savefig(out, bbox_inches="tight", facecolor=SURFACE)
    print(f"wrote {out}")


def mc2v3ScatterMatrices():
    """Total scattering matrices of both libraries, differenced, for two nuclides."""
    ref, omc, weight, upper, lower = _mc2v3Libraries()
    keep = np.flatnonzero(weight >= 1e-5)
    n = keep.size
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.6), dpi=160)
    image = None
    for ax, (label, title) in zip(axes, [("U238AA", "U-238"), ("FE56AA", "Fe-56")], strict=True):
        _style(ax, fig)
        a = np.asarray(ref[label].micros.getTotalScatterMatrix().todense())[np.ix_(keep, keep)]
        b = np.asarray(omc[label].micros.getTotalScatterMatrix().todense())[np.ix_(keep, keep)]
        # [to, from]: difference relative to the source group's total scattering in MC2-3
        scale = a.sum(axis=0)
        diff = 100.0 * (b - a) / np.where(scale > 0, scale, np.nan)[None, :]
        mask = (np.abs(a) + np.abs(b)) < 1e-6 * np.nanmax(a)
        diff[mask] = np.nan
        image = ax.imshow(diff, cmap=DIVERGING, vmin=-10, vmax=10, interpolation="nearest")
        ax.set_facecolor(SURFACE)
        ticks = list(range(0, n, 4))
        ax.set_xticks(ticks, [str(keep[t] + 1) for t in ticks])
        ax.set_yticks(ticks, [str(keep[t] + 1) for t in ticks])
        ax.set_xlabel("from group (fast → slow)", color=INK_2, fontsize=8)
        ax.set_ylabel("to group", color=INK_2, fontsize=8)
        inGroup = np.nanmean(np.abs(np.diag(diff)))
        ax.set_title(f"{title} scattering, OpenMC − MC2-3", color=INK, fontsize=9, loc="left")
        ax.text(
            0.02,
            0.03,
            f"mean |in-group diff| {inGroup:.1f}%",
            transform=ax.transAxes,
            color=INK_2,
            fontsize=7,
        )
    bar = fig.colorbar(image, ax=axes, fraction=0.025, pad=0.02, extend="both")
    bar.outline.set_visible(False)
    bar.ax.tick_params(colors=INK_2, labelsize=8)
    bar.set_label("% of the source group's total scattering", color=INK_2, fontsize=8)
    out = STATIC / "xs-mc2v3-scatter.png"
    fig.savefig(out, bbox_inches="tight", facecolor=SURFACE)
    print(f"wrote {out}")


if __name__ == "__main__":
    STATIC.mkdir(parents=True, exist_ok=True)
    c5g7PinPower()
    fftfAssemblyPower()
    mc2v3CrossSections()
    mc2v3ScatterMatrices()
