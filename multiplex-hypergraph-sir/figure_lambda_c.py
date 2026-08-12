"""Figure 2: lambda_c by subcritical-final-state extrapolation vs rho(N)=1."""
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# ---------------------------------------------------------------- style
INK, SEC, MUTED = "#1a1a1a", "#52514e", "#8a8985"
GRID, SPINE = "#e9e8e3", "#bdbcb6"
# fresh, colorblind-safe palette (validated: worst all-pairs CVD dE 11.7,
# normal-vision dE 23.4); thick lines + markers + legend supply the relief
# for the sub-3:1 contrast of teal/coral.
TEAL, CORAL, PERI = "#2BB3A3", "#F0785A", "#6C6FE0"
PT = "#4B4FA6"   # single-series accent for panels (b) and (c)

mpl.rcParams.update({
    "figure.dpi": 140, "savefig.dpi": 300,
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"], "font.size": 9,
    "mathtext.fontset": "dejavusans",
    "axes.edgecolor": SPINE, "axes.linewidth": 0.8,
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": SEC, "ytick.color": SEC,
    "xtick.labelcolor": SEC, "ytick.labelcolor": SEC,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "legend.frameon": False, "axes.axisbelow": True,
})

data = json.load(open("figure_data.json"))
by = {r["name"]: r for r in data}

fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(10.2, 3.25))
fig.subplots_adjust(left=0.065, right=0.985, bottom=0.16, top=0.90, wspace=0.34)

def panel_tag(ax, s):
    ax.text(-0.20, 1.06, s, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="bottom", ha="left", color=INK)

# ============================================================ (a) method
# plot in normalised lambda/lambda_c so the three configs share an axis and
# their linear approach to zero crosses together at the rho(N)=1 point (x=1).
show = [("2-layer m=(3,4)", TEAL,  "2-layer $m{=}(3,4)$"),
        ("6-reg. pairwise", CORAL, "6-reg. pairwise"),
        ("m=3 deg{2,3}",    PERI,  "1-layer $m{=}3$")]
# short vertical marker at lambda_c (kept low so it clears the legend above)
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
    fitx = np.linspace(xn[0], xext, 50)               # data -> zero-crossing
    axA.plot(fitx, a + b * fitx * lc_rho, "--", color=col, lw=1.0,
             alpha=0.85, zorder=2)
    fxs = np.linspace(xn[0], xn[-1], 20)
    line, = axA.plot(fxs, a + b * fxs * lc_rho, "-", color=col, lw=1.7,
                     zorder=3, marker="o", markevery=[0], ms=4.3,
                     mfc=col, mec="white", mew=0.6, label=lab)
    axA.plot(xn, ys, "o", ms=4.3, mfc=col, mec="white", mew=0.6, zorder=4)
    handles.append(line)
axA.set_xlim(0.876, 1.026)
axA.set_ylim(-0.009, 0.098)
axA.set_xlabel(r"$\lambda/\lambda_c\;[\rho(N)=1]$")
axA.set_ylabel(r"$\varepsilon/R(\infty)$  (inverse outbreak size)")
axA.xaxis.set_major_locator(MultipleLocator(0.05))
axA.text(1.0, 0.058, r"$\lambda_c$", fontsize=9, color=SEC,
         ha="center", va="bottom")
# labels listed separately, inside the panel in the empty upper-right region
# above the descending curves
axA.legend(handles=handles, loc="upper right", frameon=False, fontsize=7.5,
           handlelength=1.4, handletextpad=0.5, labelspacing=0.4,
           borderaxespad=0.6)
panel_tag(axA, "a")

# ============================================================ (b) agreement
lcr = np.array([r["lc_rho"] for r in data])
lce = np.array([r["lc_extrap"] for r in data])
sig = np.array([r["sigma"] for r in data])
lo, hi = 0.06, 0.42
axB.plot([lo, hi], [lo, hi], "-", color=MUTED, lw=1.0, zorder=1)
axB.errorbar(lcr, lce, yerr=sig, fmt="D", ms=4.6, mfc=PT, mec="white",
             mew=0.6, ecolor=PT, elinewidth=1.0, capsize=2, zorder=3)
axB.set_xlim(lo, hi); axB.set_ylim(lo, hi)
axB.set_aspect("equal")
axB.set_xlabel(r"$\lambda_c$ from $\rho(N)=1$")
axB.set_ylabel(r"$\lambda_c$ from subcritical extrap.")
axB.xaxis.set_major_locator(MultipleLocator(0.1))
axB.yaxis.set_major_locator(MultipleLocator(0.1))
maxrel = np.max(np.abs(lce - lcr) / lcr)
axB.text(0.95, 0.06, f"max rel. dev.\n{maxrel:.0e}", transform=axB.transAxes,
         fontsize=8.2, color=SEC, ha="right", va="bottom")
panel_tag(axB, "b")

# ============================================================ (c) pull / sigma
pull = (lce - lcr) / sig
order = np.argsort(lcr)
x = lcr[order]; y = pull[order]
axC.axhspan(-2, 2, color=GRID, zorder=0)
axC.axhspan(-1, 1, color="#dcdbd4", zorder=0)
axC.axhline(0, color=MUTED, lw=0.9, zorder=1)
axC.plot(x, y, "o", ms=5.0, mfc=PT, mec="white", mew=0.7, zorder=3)
axC.set_ylim(-3, 3)
axC.set_xlim(0.06, 0.42)
axC.set_xlabel(r"$\lambda_c$")
axC.set_ylabel(r"$(\lambda_c^{\mathrm{extr}}-\lambda_c^{\rho})/\sigma$")
axC.xaxis.set_major_locator(MultipleLocator(0.1))
axC.set_yticks([-2, -1, 0, 1, 2])
axC.text(0.5, 0.93, r"$\pm1\sigma,\ \pm2\sigma$", transform=axC.transAxes,
         fontsize=8.0, color=SEC, ha="center", va="top")
axC.text(0.95, 0.06, f"max\n{np.max(np.abs(pull)):.1f}$\\sigma$",
         transform=axC.transAxes, fontsize=8.5, color=SEC,
         ha="right", va="bottom")
panel_tag(axC, "c")

for ax in (axA, axB, axC):
    ax.tick_params(length=3)

fig.savefig("figure2_lambda_c.pdf", bbox_inches="tight")
fig.savefig("figure2_lambda_c.png", bbox_inches="tight", dpi=200)
print("saved figure2_lambda_c.pdf / .png")
print(f"max pull = {np.max(np.abs(pull)):.2f} sigma ; max rel dev = {maxrel:.2e}")
