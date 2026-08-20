"""Figure 3: lambda_c by subcritical-final-state extrapolation vs rho(N)=1.
(Figure 1 = mechanism schematic, Figure 2 = time evolution, both in Part 1.)

Refined journal styling: STIX serif typography, hairline chrome, minor ticks,
a colour-blind-safe jewel-tone palette (validated).
"""
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

# ---------------------------------------------------------------- palette
INK, SEC, MUTED = "#20201e", "#565550", "#9b9a93"
GRID, SPINE = "#ecebe5", "#c6c5bf"
TEAL, CORAL, INDIGO = "#1C9B8E", "#E76F51", "#4F8AC9"      # teal / coral / fresh blue
PT = "#3E6BA6"                                             # clean blue single-series accent
BAND2, BAND1 = "#f0efe9", "#e4e3db"                        # +-2sigma / +-1sigma

mpl.rcParams.update({
    "figure.dpi": 140, "savefig.dpi": 340,
    "font.family": "STIXGeneral", "mathtext.fontset": "stix",
    "font.size": 9.5,
    "axes.edgecolor": SPINE, "axes.linewidth": 0.7,
    "axes.labelcolor": INK, "axes.labelpad": 3.5, "text.color": INK,
    "xtick.color": SEC, "ytick.color": SEC,
    "xtick.labelcolor": SEC, "ytick.labelcolor": SEC,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.size": 3.4, "ytick.major.size": 3.4,
    "xtick.minor.size": 1.9, "ytick.minor.size": 1.9,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.minor.width": 0.6, "ytick.minor.width": 0.6,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.55,
    "legend.frameon": False, "axes.axisbelow": True,
    "lines.solid_capstyle": "round",
})

data = json.load(open("figure_data.json"))
by = {r["name"]: r for r in data}

fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(10.4, 3.3))
fig.subplots_adjust(left=0.062, right=0.987, bottom=0.165, top=0.905, wspace=0.335)

def tag(ax, s):
    ax.text(-0.205, 1.045, s, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="bottom", ha="left", color=INK, family="STIXGeneral")

def minorticks(ax, nx=5, ny=2):
    ax.xaxis.set_minor_locator(AutoMinorLocator(nx))
    ax.yaxis.set_minor_locator(AutoMinorLocator(ny))
    ax.tick_params(which="both", top=False, right=False)

# ============================================================ (a) method
show = [("2-layer m=(3,4)", TEAL,   "2-layer $m=(3,4)$"),
        ("6-reg. pairwise", CORAL,  "6-reg. pairwise"),
        ("m=3 deg{2,3}",    INDIGO, r"1-layer $m=3$, $k\in\{2,3\}$")]
axA.plot([1.0, 1.0], [0.0, 0.052], color=MUTED, lw=0.9, ls=(0, (4, 3)), zorder=1)
axA.axhline(0, color=MUTED, lw=0.8, zorder=1)
handles = []
for name, col, lab in show:
    r = by[name]
    lc_rho = r["lc_rho"]
    xn = np.array(r["lams"]) / lc_rho
    ys = np.array(r["ys"])
    a, b = r["fit"]
    xext = r["lc_extrap"] / lc_rho
    axA.plot(np.linspace(xn[0], xext, 60), a + b * np.linspace(xn[0], xext, 60) * lc_rho,
             "--", color=col, lw=1.0, alpha=0.75, zorder=2, dashes=(4, 2.4))
    fxs = np.linspace(xn[0], xn[-1], 40)
    line, = axA.plot(fxs, a + b * fxs * lc_rho, "-", color=col, lw=1.8, zorder=3,
                     marker="o", markevery=[0], ms=4.2, mfc=col, mec="white",
                     mew=0.7, label=lab)
    axA.plot(xn, ys, "o", ms=4.2, mfc=col, mec="white", mew=0.7, zorder=4)
    handles.append(line)
# hollow marker at the rho(N)=1 prediction, i.e. x = 1 exactly. The three fits
# cross zero at 1.000616 / 1.000386 / 1.000336 -- to the right of it and of each
# other, but by less than a pixel on this axis.
axA.plot([1.0], [0.0], "o", ms=6.5, mfc="white", mec=SEC, mew=1.0, zorder=5)
axA.set_xlim(0.876, 1.026)
axA.set_ylim(-0.009, 0.098)
axA.set_xlabel(r"$\lambda/\lambda_c\;\;[\,\rho(N)=1\,]$")
axA.set_ylabel(r"$\varepsilon/R(\infty)$   (inverse outbreak size)")
axA.xaxis.set_major_locator(MultipleLocator(0.05))
minorticks(axA)
axA.text(1.0, 0.058, r"$\lambda_c$", fontsize=10.5, color=SEC, ha="center", va="bottom")
axA.legend(handles=handles, loc="upper right", frameon=False, fontsize=8.4,
           handlelength=1.5, handletextpad=0.55, labelspacing=0.45, borderaxespad=0.7)
tag(axA, "a")

# ============================================================ (b) agreement
lcr = np.array([r["lc_rho"] for r in data])
lce = np.array([r["lc_extrap"] for r in data])
sig = np.array([r["sigma"] for r in data])
lo, hi = 0.06, 0.42
axB.plot([lo, hi], [lo, hi], "-", color=MUTED, lw=1.0, zorder=1)
axB.errorbar(lcr, lce, yerr=sig, fmt="D", ms=4.7, mfc=PT, mec="white", mew=0.7,
             ecolor=PT, elinewidth=1.0, capsize=2.2, capthick=0.8, zorder=3)
axB.set_xlim(lo, hi); axB.set_ylim(lo, hi); axB.set_aspect("equal")
axB.set_xlabel(r"$\lambda_c$   from $\rho(N)=1$")
axB.set_ylabel(r"$\lambda_c$   from subcritical extrap.")
axB.xaxis.set_major_locator(MultipleLocator(0.1))
axB.yaxis.set_major_locator(MultipleLocator(0.1))
minorticks(axB, nx=2, ny=2)
axB.text(0.145, 0.205, "identity", transform=axB.transAxes, fontsize=8.2,
         color=MUTED, ha="left", va="bottom", rotation=45,
         rotation_mode="anchor", style="italic")
maxrel = np.max(np.abs(lce - lcr) / lcr)
_e = int(np.floor(np.log10(maxrel))); _m = maxrel / 10 ** _e
if round(_m, 1) >= 10:            # keep the mantissa in [1,10)
    _m, _e = _m / 10, _e + 1
axB.text(0.965, 0.055, "max rel. dev.\n" + rf"${_m:.1f}\times10^{{{_e}}}$", transform=axB.transAxes,
         fontsize=8.6, color=SEC, ha="right", va="bottom", linespacing=1.4)
tag(axB, "b")

# ============================================================ (c) pull / sigma
pull = (lce - lcr) / sig
order = np.argsort(lcr)
x, y = lcr[order], pull[order]
axC.axhspan(-2, 2, color=BAND2, lw=0, zorder=0)
axC.axhspan(-1, 1, color=BAND1, lw=0, zorder=0)
axC.axhline(0, color=MUTED, lw=0.9, zorder=1)
axC.plot(x, y, "o", ms=5.0, mfc=PT, mec="white", mew=0.8, zorder=3)
axC.set_ylim(-3, 3); axC.set_xlim(0.06, 0.42)
axC.set_xlabel(r"$\lambda_c$")
axC.set_ylabel(r"$(\lambda_c^{\mathrm{extr}}-\lambda_c^{\rho})\,/\,\sigma$")
axC.xaxis.set_major_locator(MultipleLocator(0.1))
axC.set_yticks([-2, -1, 0, 1, 2])
minorticks(axC, nx=2, ny=1)
axC.text(0.5, 0.955, r"$\pm1\sigma,\;\;\pm2\sigma$", transform=axC.transAxes,
         fontsize=8.6, color=SEC, ha="center", va="top")
axC.text(0.965, 0.055, "max\n" + rf"${np.max(np.abs(pull)):.2f}\,\sigma$", transform=axC.transAxes,
         fontsize=8.8, color=SEC, ha="right", va="bottom", linespacing=1.4)
tag(axC, "c")

fig.savefig("figure3_lambda_c.pdf", bbox_inches="tight")
fig.savefig("figure3_lambda_c.png", bbox_inches="tight", dpi=210)
print("saved; max pull = %.2f sigma ; max rel dev = %.2e"
      % (np.max(np.abs(pull)), maxrel))
