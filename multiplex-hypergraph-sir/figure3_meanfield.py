"""
Figure 3: turning off the three mean-field deviations term by term, and the
sign and magnitude of the net deviation on the (m, lambda) plane.
"""
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

from theory import lambda_c_rho
from meanfield import factors, lambda_c_mf, lambda_c_switched

# ---------------------------------------------------------------- style
INK, SEC, MUTED = "#20201e", "#565550", "#9b9a93"
GRID, SPINE = "#ecebe5", "#c6c5bf"
TEAL, CORAL, INDIGO = "#2BB3A3", "#F0785A", "#6C6FE0"
PT = "#5960CE"
mpl.rcParams.update({
    "figure.dpi": 140, "savefig.dpi": 340,
    "font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 9.5,
    "axes.edgecolor": SPINE, "axes.linewidth": 0.7,
    "axes.labelcolor": INK, "axes.labelpad": 3.5, "text.color": INK,
    "xtick.color": SEC, "ytick.color": SEC,
    "xtick.labelcolor": SEC, "ytick.labelcolor": SEC,
    "xtick.major.size": 3.4, "ytick.major.size": 3.4,
    "xtick.minor.size": 1.9, "ytick.minor.size": 1.9,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.55,
    "legend.frameon": False, "axes.axisbelow": True,
    "lines.solid_capstyle": "round",
})
# diverging: blue <-> neutral grey <-> red  (never a hue at the midpoint)
DIV = LinearSegmentedColormap.from_list(
    "div", ["#0E7F72", "#3FB3A4", "#A6DBD2", "#F2F0EA", "#FBD2C1", "#F08D6C", "#C2542F"])

P = {(3, 3): 0.5, (5, 5): 0.5}          # k in {3,5}, two layers, m1 = m2 = m
MS = list(range(2, 17))
LAMS = np.geomspace(0.01, 2.0, 260)

# ---------------------------------------------------------------- data
Z = np.empty((len(MS), len(LAMS)))       # log2(rho_MF / rho_exact)
for i, m in enumerate(MS):
    for j, lam in enumerate(LAMS):
        Z[i, j] = -np.log2(factors(P, (m, m), lam)[3])   # net = rho/rho_MF

lc_ex = np.array([lambda_c_rho(P, (m, m)) for m in MS])
lc_mf = np.array([lambda_c_mf(P, (m, m)) for m in MS])
ratio = lc_mf / lc_ex

M_A = 10                                  # representative m for panel (a)
fT = np.array([factors(P, (M_A, M_A), l)[0] for l in LAMS])
fD = np.array([factors(P, (M_A, M_A), l)[1] for l in LAMS])
fC = np.array([factors(P, (M_A, M_A), l)[2] for l in LAMS])
net = fT * fD * fC

# ---------------------------------------------------------------- figure
fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(10.4, 3.3))
fig.subplots_adjust(left=0.062, right=0.975, bottom=0.165, top=0.905, wspace=0.40)

def tag(ax, s):
    ax.text(-0.205, 1.045, s, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="bottom", ha="left", color=INK)

# ---- (a) the three factors, switched on one at a time -------------------
axA.axhline(1.0, color=MUTED, lw=0.9, ls=(0, (4, 3)), zorder=1)
axA.plot(LAMS, fC, "-", color=TEAL, lw=1.8, zorder=3, label=r"$f_C$  cascade")
axA.plot(LAMS, fD * np.ones_like(LAMS), "-", color=INDIGO, lw=1.8, zorder=3,
         label=r"$f_D$  excess $-\delta_{ab}$")
axA.plot(LAMS, fT, "-", color=CORAL, lw=1.8, zorder=3, label=r"$f_T$  $\lambda\!\to\!T$")
axA.plot(LAMS, net, "-", color=INK, lw=2.1, zorder=4, label=r"net  $=f_Tf_Df_C$")
above = net > 1
if above.any():
    axA.fill_between(LAMS, 1, net, where=above, color=TEAL, alpha=0.16,
                     lw=0, zorder=2)
axA.set_xscale("log")
axA.set_xlim(LAMS[0], LAMS[-1]); axA.set_ylim(0.25, 2.02)
axA.set_xlabel(r"$\lambda$")
axA.set_ylabel(r"factor on $\rho$   (exact / mean field)")
axA.yaxis.set_minor_locator(AutoMinorLocator(2))
axA.tick_params(which="both", top=False, right=False)
axA.legend(loc="lower left", frameon=False, fontsize=7.6, handlelength=1.4,
           handletextpad=0.5, labelspacing=0.35, borderaxespad=0.55)
axA.text(0.97, 0.955, rf"$m={M_A}$", transform=axA.transAxes, fontsize=9,
         color=SEC, ha="right", va="top")
tag(axA, "a")

# ---- (b) net deviation on the (m, lambda) plane -------------------------
vmax = float(np.max(np.abs(Z)))
norm = TwoSlopeNorm(vmin=-max(0.75, np.max(-Z)), vcenter=0.0, vmax=vmax)
pcm = axB.pcolormesh(LAMS, MS, Z, cmap=DIV, norm=norm, shading="nearest")
cs = axB.contour(LAMS, MS, Z, levels=[0.0], colors=[INK], linewidths=1.4,
                 linestyles="-")
axB.clabel(cs, fmt={0.0: "0"}, fontsize=8, inline=True)
axB.plot(lc_ex, MS, "-", color="white", lw=2.6, zorder=4)
axB.plot(lc_ex, MS, "--", color=INK, lw=1.3, zorder=5)
axB.text(lc_ex[2] * 1.55, MS[2], r"$\lambda_c$", color=INK, fontsize=9.5,
         ha="left", va="center", zorder=6)
axB.set_xscale("log")
axB.set_xlim(LAMS[0], LAMS[-1]); axB.set_ylim(MS[0] - 0.5, MS[-1] + 0.5)
axB.set_xlabel(r"$\lambda$")
axB.set_ylabel(r"group size $m$")
axB.yaxis.set_major_locator(MultipleLocator(2))
axB.grid(False)
axB.tick_params(which="both", top=False, right=False)
cb = fig.colorbar(pcm, ax=axB, pad=0.03, aspect=22)
cb.set_label(r"$\log_2(\rho^{\mathrm{MF}}/\rho)$", fontsize=8.6)
cb.ax.tick_params(labelsize=7.8, color=SEC, labelcolor=SEC)
cb.outline.set_edgecolor(SPINE); cb.outline.set_linewidth(0.6)
axB.text(0.60, 0.90, "MF under-\nestimates $\\rho$", transform=axB.transAxes,
         fontsize=7.8, color="#0A5D53", ha="center", va="top", linespacing=1.35)
axB.text(0.955, 0.10, "MF over-\nestimates $\\rho$", transform=axB.transAxes,
         fontsize=7.8, color="#8A3A1E", ha="right", va="bottom", linespacing=1.35)
tag(axB, "b")

# ---- (c) the threshold consequence --------------------------------------
axC.axhline(1.0, color=MUTED, lw=0.9, ls=(0, (4, 3)), zorder=1)
axC.plot(MS, ratio, "-o", color=PT, lw=1.7, ms=4.4, mfc=PT, mec="white",
         mew=0.7, zorder=3)
axC.set_xlim(MS[0] - 0.4, MS[-1] + 0.4); axC.set_ylim(0.72, 1.03)
axC.set_xlabel(r"group size $m$")
axC.set_ylabel(r"$\lambda_c^{\mathrm{MF}}/\lambda_c$")
axC.xaxis.set_major_locator(MultipleLocator(2))
axC.yaxis.set_major_locator(MultipleLocator(0.05))
axC.xaxis.set_minor_locator(AutoMinorLocator(2))
axC.yaxis.set_minor_locator(AutoMinorLocator(2))
axC.tick_params(which="both", top=False, right=False)
for m, dx, dy in [(2, 8, -11), (16, -2, -13)]:
    r = ratio[MS.index(m)]
    axC.annotate(f"{r:.2f}", xy=(m, r), xytext=(dx, dy),
                 textcoords="offset points", fontsize=7.8, color=SEC,
                 ha="center")
axC.text(0.955, 0.26, "mean field always\nunderestimates $\\lambda_c$",
         transform=axC.transAxes, fontsize=8.2, color=SEC, ha="right",
         va="center", linespacing=1.4)
tag(axC, "c")

fig.savefig("figure3_meanfield.pdf", bbox_inches="tight")
fig.savefig("figure3_meanfield.png", bbox_inches="tight", dpi=210)

# ---------------------------------------------------------------- data out
flip = {}
for m in MS:
    ls = np.geomspace(1e-3, 20, 3000)
    nets = np.array([factors(P, (m, m), l)[3] for l in ls])
    a = ls[nets > 1]
    flip[m] = [float(a.min()), float(a.max())] if len(a) else None
json.dump({"m": MS, "lc_exact": lc_ex.tolist(), "lc_mf": lc_mf.tolist(),
           "ratio": ratio.tolist(), "flip_window": flip},
          open("figure3_data.json", "w"), indent=2)

print("saved figure3_meanfield.pdf/.png")
print("lam_c^MF/lam_c :", ", ".join(f"m={m}:{r:.4f}" for m, r in zip(MS, ratio) if m in (2,3,4,6,8,16)))
print("first m with sign flip:", min(m for m in MS if flip[m]))
print("flip window at m=8 :", flip[8])
print("max log2(rho_MF/rho):", f"{Z.max():.3f}", " min:", f"{Z.min():.3f}")
