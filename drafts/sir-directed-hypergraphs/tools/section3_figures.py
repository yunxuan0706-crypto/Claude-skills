#!/usr/bin/env python3
"""Section 3 figures, to SIAM Review specification.

Layout follows the target journal rather than a generic preset. SIAM Review sets
a single 5.57 in text column (measured from the journal PDF: 486 x 720 pt page,
text block 58.9-459.7 pt), body face Times, and figure text must stay legible at
8 pt in print. Figures are drawn at 5.5 in so they sit at 100% scale in that
column -- a figure that has to be shrunk to fit loses its type size with it.

Two figures rather than one four-panel block, because they make two different
arguments:

  Fig. 1  the closure reproduces the whole time course, not just the endpoints;
  Fig. 2  and it is quantitatively right about finite size, threshold location,
          and structural dependence.

Colours are Okabe-Ito, with line style and marker carrying the same information
redundantly so the panels survive greyscale printing and colour-vision
deficiency.
"""
import json
import math
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SKILL = "/root/.claude/skills/scipilot-figure-skill/scripts"
sys.path.insert(0, SKILL)
from export_figure import export_figure          # noqa: E402
from layout_tools import add_panel_labels, finalize_figure   # noqa: E402
import visual_qa                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "section3.json")
FIGS = os.path.join(HERE, "..", "figures")

# Okabe-Ito, colour-vision-deficiency safe
C_SIM = "#000000"     # simulation: black
C_EBCM = "#0072B2"    # blue
C_MF = "#D55E00"      # vermillion
C_REF = "#999999"

def eq_number(label, src=os.path.join(HERE, "..", "section3.src.md"), sec=3):
    """Resolve \\label{...} to the number the manuscript build gives it.

    A hard-coded "Eq. (38)" in a legend goes stale the moment an equation is
    inserted upstream -- and unlike a \\eqref in the text, nothing checks it.
    """
    n = 0
    for m in re.finditer(r"^\$\$\n(.*?)^\$\$", open(src, encoding="utf-8").read(),
                         re.S | re.M):
        if "\\notag" in m.group(1):
            continue
        n += 1
        lm = re.search(r"\\label\{([^}]+)\}", m.group(1))
        if lm and lm.group(1) == label:
            return f"{sec}.{n}"
    raise SystemExit(f"figure legend references an unknown label: {label}")


WIDTH = 5.5          # SIAM Review text column is 5.57 in

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


def use_cjk(on):
    """Chinese labels for the draft; English for submission."""
    if on:
        plt.rcParams["font.family"] = "serif"
        plt.rcParams["font.serif"] = ["Noto Serif CJK SC", "Liberation Serif"]
    else:
        plt.rcParams["font.serif"] = ["Liberation Serif", "DejaVu Serif"]


def L(zh, en, cn):
    return zh if cn else en



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


# --------------------------------------------------------------- figure 1 ----
def figure1(d, cn):
    tr = d["trajectory"]
    t = np.array(tr["t"])
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 2.1))

    for ax, key, ylab in (
        (axes[0], "S", L("易感比例 $S(t)$", "Susceptible fraction $S(t)$", cn)),
        (axes[1], "I", L("感染比例 $I(t)$", "Infected fraction $I(t)$", cn)),
    ):
        sim = np.array(tr[f"{key}_sim"])
        sem = np.array(tr[f"{key}_sem"])
        ax.fill_between(t, sim - 1.96 * sem, sim + 1.96 * sem,
                        color=C_SIM, alpha=0.22, lw=0, zorder=2)
        ax.plot(t, sim, color=C_SIM, lw=1.1, zorder=3,
                label=L("精确仿真", "Exact simulation", cn))
        ax.plot(t, tr[f"{key}_ebcm"], color=C_EBCM, ls="--", lw=1.1, zorder=4,
                label=L("并发型 EBCM", "Concurrent EBCM", cn))
        ax.plot(t, tr[f"{key}_mf"], color=C_MF, ls=":", lw=1.2, zorder=4,
                label=L("均场闭合", "Mean field", cn))
        ax.set_xlabel(L("时间 $t$（以 $1/\\mu$ 为单位）",
                        "Time $t$ (units of $1/\\mu$)", cn))
        ax.set_ylabel(ylab)
        ax.set_xlim(0, t.max())
        ax.margins(y=0.06)

    axes[1].set_ylim(bottom=0)
    axes[0].legend(loc="lower left", handlelength=2.2)
    txt = L(f"$N={tr['N']}$, $\\tau=\\eta=2$, $\\lambda={tr['ratio']}\\lambda_c$",
            f"$N={tr['N']}$, $\\tau=\\eta=2$, $\\lambda={tr['ratio']}\\lambda_c$", cn)
    axes[1].text(0.97, 0.93, txt, transform=axes[1].transAxes,
                 ha="right", va="top", fontsize=7)
    return fig


# --------------------------------------------------------------- figure 2 ----
def figure2(d, cn):
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 2.45))

    # (a) closure residual versus N
    ax = axes[0]
    fs = d["finite_size"]
    N = np.array([r["N"] for r in fs], float)
    for key, col, mk, ls, lab in (
        ("ebcm", C_EBCM, "o", "-", L("并发型 EBCM", "Concurrent EBCM", cn)),
        ("mf", C_MF, "s", "--", L("均场闭合", "Mean field", cn)),
    ):
        y = np.array([r[key] for r in fs])
        e = np.array([r[key + "_sem"] for r in fs])
        ax.errorbar(N, y, yerr=e, color=col, marker=mk, ls=ls, capsize=1.8,
                    elinewidth=0.6, capthick=0.6, label=lab)
    # Fit the EBCM decay rather than drawing an assumed N^-1 guide: the data
    # are steeper than N^-1, so the guide would have misrepresented them.
    ye = np.array([r["ebcm"] for r in fs])
    ee = np.array([r["ebcm_sem"] for r in fs])
    w = 1.0 / (ee / ye) ** 2                       # weights on log residuals
    lx, ly = np.log(N), np.log(ye)
    Sw, Sx, Sy = w.sum(), (w * lx).sum(), (w * ly).sum()
    Sxx, Sxy = (w * lx * lx).sum(), (w * lx * ly).sum()
    det = Sw * Sxx - Sx * Sx
    slope = (Sw * Sxy - Sx * Sy) / det
    icpt = (Sxx * Sy - Sx * Sxy) / det
    dslope = math.sqrt(Sw / det)
    xs = np.logspace(math.log10(N[0]), math.log10(N[-1]), 50)
    ax.plot(xs, np.exp(icpt) * xs ** slope, color=C_REF, ls="-.", lw=0.7,
            zorder=1,
            label=L(f"$\\propto N^{{{slope:.1f}}}$", f"$\\propto N^{{{slope:.1f}}}$", cn))
    print(f"  [fit] EBCM residual ~ N^({slope:.3f} +/- {dslope:.3f})")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks(N)
    ax.set_xticklabels([f"{int(v)}" for v in N])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlabel(L("系统规模 $N$", "System size $N$", cn))
    ax.set_ylabel(L("$S(t)$ 平均绝对偏差", "Mean absolute error in $S(t)$", cn))
    ax.legend(loc="lower left", handlelength=2.0)

    # (b) final size across the threshold
    ax = axes[1]
    th = d["threshold"]
    r = np.array([x["ratio"] for x in th["rows"]])
    sim = np.array([x["sim"] for x in th["rows"]])
    sem = np.array([x["sim_sem"] for x in th["rows"]])
    ax.axvline(1.0, color=C_REF, lw=0.6, ls="-")
    ax.axvline(th["lc_mf"] / th["lc"], color=C_MF, lw=0.6, ls=":")
    ax.errorbar(r, sim, yerr=1.96 * sem, color=C_SIM, marker="o", ls="none",
                capsize=1.8, elinewidth=0.6, capthick=0.6, zorder=5,
                label=L("精确仿真", "Exact simulation", cn))
    ax.plot(r, [x["ebcm"] for x in th["rows"]], color=C_EBCM, ls="--",
            marker="", lw=1.1, label=L("并发型 EBCM", "Concurrent EBCM", cn))
    ax.plot(r, [x["mf"] for x in th["rows"]], color=C_MF, ls=":", lw=1.2,
            label=L("均场闭合", "Mean field", cn))
    ax.set_xlabel(L("$\\lambda/\\lambda_c$", "$\\lambda/\\lambda_c$", cn))
    ax.set_ylabel(L("终态规模 $R(\\infty)$", "Final size $R(\\infty)$", cn))
    ax.set_ylim(0, 0.98)
    # Name the two verticals in the legend. Placing the labels in the axes ran
    # them into the y-axis on the left and into the data on the right. In a
    # 1.8 in panel a four-row legend spans the full width whatever its corner,
    # so the y-limit reserves a blank band above the curves for it to sit in.
    from matplotlib.lines import Line2D
    handles, labels = ax.get_legend_handles_labels()
    handles += [Line2D([], [], color=C_MF, ls=":", lw=0.8),
                Line2D([], [], color=C_REF, ls="-", lw=0.8)]
    labels += [L("均场阈值 $\\lambda_c^{\\rm MF}$", "MF threshold $\\lambda_c^{\\rm MF}$", cn),
               L("真实阈值 $\\lambda_c$", "True threshold $\\lambda_c$", cn)]
    ax.legend(handles, labels, loc="upper left", handlelength=1.9,
              labelspacing=0.32, borderpad=0.3)

    # (c) threshold versus in-out degree correlation
    LCRIO = eq_number("eq:lcrio")
    ax = axes[2]
    ri = d["rio"]
    x = np.array([q["rio"] for q in ri])
    y = np.array([q["lc"] for q in ri])
    xe = np.array([q["rio_sem"] for q in ri])
    ye = np.array([q["lc_sem"] for q in ri])
    xs = np.linspace(x.min(), x.max(), 200)
    ki = np.mean([q["ki_mean"] for q in ri])
    si = np.mean([q["si"] for q in ri])
    so = np.mean([q["so"] for q in ri])
    nm = np.mean([q["NoverM"] for q in ri])
    ax.plot(xs, 1.0 / (2 * ki + nm * xs * si * so - 1), color=C_EBCM, lw=1.0,
            label=L(f"式 ({LCRIO})", f"Eq. ({LCRIO})", cn))
    ax.errorbar(x, y, xerr=xe, yerr=ye, color=C_SIM, marker="D", ls="none",
                capsize=1.8, elinewidth=0.6, capthick=0.6, zorder=5,
                label=L("位形模型实现", "Realisations", cn))
    ax.set_xlabel(L("入出度相关 $r_{io}$",
                    "In-out degree correlation $r_{io}$", cn))
    ax.set_ylabel(L("爆发阈值 $\\lambda_c$", "Epidemic threshold $\\lambda_c$", cn))
    ax.legend(loc="upper right", handlelength=1.8)
    return fig


def main(cn=True, tag="zh"):
    d = json.load(open(DATA))
    os.makedirs(FIGS, exist_ok=True)
    plt.rcParams.update(SIREV)
    use_cjk(cn)

    for name, builder in (("fig1_trajectory", figure1), ("fig2_validation", figure2)):
        fig = builder(d, cn)
        finalize_figure(fig)
        check_legends(fig)
        add_panel_labels(fig, style="nature")
        base = os.path.join(FIGS, f"{name}_{tag}")
        prev = visual_qa.render_preview(fig, base + "_preview.png", dpi=220)
        issues = visual_qa.audit_layout(fig)
        print(f"\n{name} [{tag}]")
        print(visual_qa.print_report(issues))
        export_figure(fig, base, formats=["pdf", "png"], dpi=600,
                      grayscale_preview=True)
        plt.close(fig)
        print("  preview:", prev)


if __name__ == "__main__":
    # One set of files, English-labelled: journal figures carry English labels
    # even when the manuscript is drafted in Chinese. main(cn=True, tag="zh")
    # regenerates the Chinese-labelled variant; nothing but the labels changes.
    main(cn=False, tag="en")
