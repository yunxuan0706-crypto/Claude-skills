#!/usr/bin/env python3
"""Figures for the verification of section 3.1, drawn from data/verification.json.

Labels are English throughout. Three guards run on every figure before it is
written: a legend must stay inside its axes, no two pieces of text may overlap,
and no text or legend may sit on top of a plotted curve. Each of the three was
added after a defect that the programmatic layout audit had passed.
"""
import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "verification.json")
FIGS = os.path.join(HERE, "..", "figures")

WIDTH = 5.5          # single text column

C_SIM = "#000000"     # exact simulation
C_CL = "#0072B2"      # group closure
C_MF = "#D55E00"      # mean field
C_REF = "#999999"

SIREV = {
    "font.family": "serif",
    "font.serif": ["Liberation Serif", "DejaVu Serif"],
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "grid.linewidth": 0.4,
    "lines.linewidth": 1.0,
    "lines.markersize": 3.2,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "legend.frameon": False,
    # Body face is Times-metric, so the maths must be too -- matplotlib's
    # default mathtext is DejaVu Sans and reads as a different document.
    "mathtext.fontset": "stix",
    "axes.unicode_minus": False,
    "savefig.bbox": "standard",
    "figure.constrained_layout.use": True,
}



def check_legends(fig):
    """Fail if a legend spills outside its axes.

    The layout audit checks glyphs, clipping and tick overlap; it does not
    notice a legend that has outgrown its panel, which is how a long label
    silently covers the data.
    """
    fig.canvas.draw()
    bad = []
    for ax in fig.axes:
        lg = ax.get_legend()
        if lg is None:
            continue
        lb = lg.get_window_extent()
        ab = ax.get_window_extent()
        if (lb.x0 < ab.x0 - 1 or lb.x1 > ab.x1 + 1
                or lb.y0 < ab.y0 - 1 or lb.y1 > ab.y1 + 1):
            bad.append(ax.get_ylabel() or ax.get_xlabel() or "axes")
    if bad:
        raise SystemExit("legend overflows its axes in: " + "; ".join(bad))


def check_text_overlap(fig, pad=0.5):
    """Fail if any two pieces of text in the figure overlap.

    Covers tick labels, axis labels, in-axes annotations, panel letters and
    legend entries, across panels as well as within them -- the cases a
    per-axes tick check cannot see. Run it after every artist is placed.
    """
    fig.canvas.draw()
    items = []                       # (label, bbox, legend-id or None)
    for i, ax in enumerate(fig.axes):
        group = [(f"panel{i} xlabel", ax.xaxis.label), (f"panel{i} ylabel", ax.yaxis.label)]
        # Matplotlib keeps tick labels for ticks outside the view; those are
        # never drawn, so comparing their boxes invents collisions.
        for axis, kind in ((ax.xaxis, "x"), (ax.yaxis, "y")):
            lo, hi = sorted(ax.get_xlim() if kind == "x" else ax.get_ylim())
            locs = axis.get_ticklocs()
            labs = axis.get_ticklabels()
            group += [(f"panel{i} {kind}tick {t.get_text()}", t)
                      for v, t in zip(locs, labs) if lo <= v <= hi]
        group += [(f"panel{i} text {t.get_text()[:18]}", t) for t in ax.texts]
        for name, t in group:
            if t.get_text().strip() and t.get_visible():
                items.append((name, t.get_window_extent(), None))
        lg = ax.get_legend()
        if lg is not None:
            for t in lg.get_texts():
                items.append((f"panel{i} legend {t.get_text()[:18]}",
                              t.get_window_extent(), id(lg)))
    for t in fig.texts:
        if t.get_text().strip() and t.get_visible():
            items.append((f"figure text {t.get_text()[:18]}", t.get_window_extent(), None))

    clashes = []
    for i in range(len(items)):
        ni, bi, gi = items[i]
        for j in range(i + 1, len(items)):
            nj, bj, gj = items[j]
            if gi is not None and gi == gj:
                continue            # rows of one legend are laid out, not colliding
            if (bi.x0 < bj.x1 - pad and bj.x0 < bi.x1 - pad
                    and bi.y0 < bj.y1 - pad and bj.y0 < bi.y1 - pad):
                clashes.append(f"{ni} <-> {nj}")
    if clashes:
        raise SystemExit("overlapping labels:\n  " + "\n  ".join(clashes))


def check_text_over_data(fig, pad=1.0, samples=24):
    """Fail if a label or legend sits on top of a plotted curve.

    check_text_overlap only compares text with text. English labels are far
    wider than their Chinese counterparts, so a legend that cleared the curves
    in one language can land squarely on them in the other -- which is exactly
    what happened to figure 1.
    """
    import numpy as _np
    fig.canvas.draw()
    boxes = []
    for i, ax in enumerate(fig.axes):
        for t in ax.texts:
            if t.get_text().strip() and t.get_visible():
                boxes.append((ax, f"panel{i} text {t.get_text()[:18]}",
                              t.get_window_extent()))
        lg = ax.get_legend()
        if lg is not None:
            boxes.append((ax, f"panel{i} legend", lg.get_window_extent()))
    hits = []
    for ax, name, bb in boxes:
        for ln in ax.lines:
            xy = ln.get_xydata()
            if len(xy) < 2:
                continue                      # axvline/axhline: guides, not data
            pix = ax.transData.transform(xy)
            # sample along each segment: a curve can cross a box between vertices
            a, b = pix[:-1], pix[1:]
            ts = _np.linspace(0, 1, samples)[:, None, None]
            pts = (a[None] + ts * (b - a)[None]).reshape(-1, 2)
            inside = ((pts[:, 0] > bb.x0 + pad) & (pts[:, 0] < bb.x1 - pad)
                      & (pts[:, 1] > bb.y0 + pad) & (pts[:, 1] < bb.y1 - pad))
            if inside.any():
                hits.append(f"{name} <-> {ln.get_label() or 'curve'}")
    if hits:
        raise SystemExit("text over data:\n  " + "\n  ".join(sorted(set(hits))))




# ------------------------------------------------------------- figure 1 ----
def figure1(d):
    """Time course: exact simulation against the closure and the mean field."""
    tr = d["trajectory"]
    t = np.array(tr["t"])
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 2.2))
    for ax, key, ylab in ((axes[0], "S", "Susceptible fraction $S(t)$"),
                          (axes[1], "I", "Infected fraction $I(t)$")):
        sim = np.array(tr[f"{key}_sim"]); sem = np.array(tr[f"{key}_sem"])
        ax.fill_between(t, sim-1.96*sem, sim+1.96*sem, color=C_SIM, alpha=0.22,
                        lw=0, zorder=2)
        ax.plot(t, sim, color=C_SIM, lw=1.1, zorder=3, label="Exact simulation")
        ax.plot(t, tr[f"{key}_closure"], color=C_CL, ls="--", lw=1.1, zorder=4,
                label="Group closure")
        ax.plot(t, tr[f"{key}_mf"], color=C_MF, ls=":", lw=1.3, zorder=4,
                label="Mean field")
        ax.set_xlabel("Time $t$ (units of $1/\\mu$)")
        ax.set_ylabel(ylab)
        ax.set_xlim(0, t.max())
        ax.margins(y=0.06)
    axes[1].set_ylim(bottom=0)
    # Both panels fall away from the upper right by t ~ 7, so that corner stays
    # clear of every curve whatever the label widths.
    axes[0].legend(loc="upper right", handlelength=2.0, labelspacing=0.3)
    s = d["setup"]
    axes[1].text(0.97, 0.95,
                 f"$N={tr['N']}$, $m=({s['ms'][0]},{s['ms'][1]})$, "
                 f"$\\lambda={s['ratio']}\\lambda_c$",
                 transform=axes[1].transAxes, ha="right", va="top", fontsize=7)
    return fig


# ------------------------------------------------------------- figure 2 ----
def figure2(d):
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 2.05))

    # (a) how the two closures' errors scale with N
    ax = axes[0]
    sc = d["scaling"]
    N = np.array([r["N"] for r in sc], float)
    for key, col, mk, ls, lab in (("closure", C_CL, "o", "--", "Group closure"),
                                  ("meanfield", C_MF, "s", ":", "Mean field")):
        ax.plot(N, [r[key] for r in sc], color=col, marker=mk, ls=ls, label=lab)
    ax.plot(N, [r["noise"] for r in sc], color=C_REF, ls="-.", lw=0.7,
            label="Sampling floor")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks(N); ax.set_xticklabels([f"{int(v)}" for v in N])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlabel("System size $N$")
    ax.set_ylabel("Mean $|\\Delta S|$")
    # Reserve a band below the smallest datum; on a log axis this costs little
    # and it is the only way a three-row legend fits without covering the data.
    lo = min(min(r["closure"] for r in sc), min(r["noise"] for r in sc))
    ax.set_ylim(bottom=lo/40)
    ax.legend(loc="lower left", handlelength=1.7, labelspacing=0.25,
              borderpad=0.25, fontsize=6)

    # (b) the within-group cascade against its Monte-Carlo check
    ax = axes[1]
    ca = d["cascade"]
    lams = sorted({r["lam"] for r in ca})
    for lam, col in zip(lams, (C_CL, C_MF, C_SIM)):
        rs = sorted([r for r in ca if r["lam"] == lam], key=lambda r: r["m"])
        if len(rs) < 2: continue
        m = [r["m"] for r in rs]
        ax.plot(m, [r["recursion"] for r in rs], color=col, lw=1.0)
        ax.errorbar(m, [r["mc"] for r in rs], yerr=[2*r["mc_se"] for r in rs],
                    color=col, marker="o", ls="none", ms=3, capsize=1.6,
                    elinewidth=0.6)
        ax.plot(m, [r["naive"] for r in rs], color=col, ls=":", lw=0.8)
    ax.set_xlabel("Group size $m$")
    ax.set_ylabel("Cascade size $C$")
    ax.set_xticks([2, 3, 4, 5, 6])
    from matplotlib.lines import Line2D
    ax.legend([Line2D([], [], color=C_REF, lw=1.0),
               Line2D([], [], color=C_REF, marker="o", ls="none", ms=3),
               Line2D([], [], color=C_REF, ls=":", lw=0.8)],
              ["Recursion", "Monte Carlo", "$(m-1)T$"],
              loc="upper left", handlelength=1.5, labelspacing=0.25,
              borderpad=0.25, fontsize=6)

    # (c) two independent routes to the threshold
    ax = axes[2]
    tc = d["threshold_crosscheck"]
    x = [r["spectral"] for r in tc]; y = [r["jacobian"] for r in tc]
    lo, hi = min(x+y)*0.6, max(x+y)*1.6
    ax.plot([lo, hi], [lo, hi], color=C_REF, ls="-", lw=0.7, zorder=1)
    ax.plot(x, y, color=C_SIM, marker="D", ls="none", ms=3.4, zorder=3,
            label="6 configurations")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("$\\lambda_c$ from $\\rho(\\mathsf{N})=1$")
    ax.set_ylabel("$\\lambda_c$ from Jacobian")
    worst = max(r["rel"] for r in tc)
    ax.text(0.95, 0.06, f"max rel. dev.\n{worst:.0e}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=6)
    return fig


def main():
    d = json.load(open(DATA))
    os.makedirs(FIGS, exist_ok=True)
    plt.rcParams.update(SIREV)
    for name, builder in (("fig1_trajectory", figure1), ("fig2_verification", figure2)):
        fig = builder(d)
        fig.tight_layout()
        for ax, lab in zip(fig.axes, "abcdefgh"):
            ax.text(-0.02, 1.06, lab, transform=ax.transAxes, fontweight="bold",
                    fontsize=9, ha="right", va="bottom")
        check_legends(fig)
        check_text_overlap(fig)
        check_text_over_data(fig)
        base = os.path.join(FIGS, name)
        for ext in ("pdf", "png"):
            fig.savefig(f"{base}.{ext}", dpi=600, bbox_inches="tight")
        print("wrote", base + ".{pdf,png}")
        plt.close(fig)


if __name__ == "__main__":
    main()
