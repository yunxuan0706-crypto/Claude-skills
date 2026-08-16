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
# chosen from the author's palette; the one factor that pushes rho UP is warm,
# the two that push it DOWN are cool (validated all-pairs: CVD dE 15.8,
# normal-vision 21.3)
C_UP, C_DN1, C_DN2 = "#FD7E21", "#1A9CFC", "#2F3DFA"
PT = "#2F3DFA"          # single series, panel (c)
SIMC = "#2F3DFA"        # simulation, panel (d)
MFC = "#FF5359"         # mean field, panel (d)
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
# author-specified scale for panel (b): blue -> light blue -> yellow -> orange
# -> red. The five stops are evenly spaced and TwoSlopeNorm puts vcenter=0 at
# the middle stop, so the yellow band sits exactly on the sign-flip boundary:
# cool = mean field underestimates rho, warm = it overestimates.
DIV = LinearSegmentedColormap.from_list(
    "div", ["#5271AE", "#70ACDE", "#F5CC7D", "#FFA660", "#D85B59"])

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
fig, (axA, axB, axC, axD) = plt.subplots(1, 4, figsize=(13.4, 3.25))
fig.subplots_adjust(left=0.049, right=0.988, bottom=0.165, top=0.905, wspace=0.42)

def tag(ax, s):
    ax.text(-0.24, 1.045, s, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="bottom", ha="left", color=INK)

# ---- (a) the three factors, switched on one at a time -------------------
axA.axhline(1.0, color=MUTED, lw=0.9, ls=(0, (4, 3)), zorder=1)
axA.plot(LAMS, fC, "-", color=C_UP, lw=1.8, zorder=3, label=r"$f_C$  cascade")
axA.plot(LAMS, fD * np.ones_like(LAMS), "-", color=C_DN1, lw=1.8, zorder=3,
         label=r"$f_D$  excess $-\delta_{ab}$")
axA.plot(LAMS, fT, "-", color=C_DN2, lw=1.8, zorder=3, label=r"$f_T$  $\lambda\!\to\!T$")
axA.plot(LAMS, net, "-", color=INK, lw=2.1, zorder=4, label=r"net  $=f_Tf_Df_C$")
above = net > 1
if above.any():
    axA.fill_between(LAMS, 1, net, where=above, color=C_UP, alpha=0.15,
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
         fontsize=7.8, color="#2F4A7A", ha="center", va="top", linespacing=1.35)
axB.text(0.955, 0.10, "MF over-\nestimates $\\rho$", transform=axB.transAxes,
         fontsize=7.8, color="#8E3230", ha="right", va="bottom", linespacing=1.35)
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

# ---- (d) simulation validation -----------------------------------------
sim = json.load(open("figure3_sim.json"))
axD.axhline(1.0, color=MUTED, lw=0.9, ls=(0, (4, 3)), zorder=1)
for frac, mk, fc in ((0.6, "o", "white"), (0.8, "s", SIMC)):
    rs = [r for r in sim if abs(r["frac"] - frac) < 1e-9]
    xs = [r["m"] for r in rs]
    ys = [r["chi_sim"] / r["chi_exact"] for r in rs]
    es = [r["sem"] / r["chi_exact"] for r in rs]
    axD.errorbar(xs, ys, yerr=es, fmt=mk, ms=5.0, mfc=fc, mec=SIMC, mew=1.2,
                 ecolor=SIMC, elinewidth=1.0, capsize=2.4, capthick=0.8,
                 zorder=4, ls="none", label=rf"sim, $\lambda={frac}\lambda_c$")
    axD.plot(xs, [r["chi_mf"] / r["chi_exact"] for r in rs], "-", color=MFC,
             lw=1.6, marker=mk, ms=4.0, mfc="white", mec=MFC, mew=1.1,
             zorder=3, label=rf"mean field, $\lambda={frac}\lambda_c$")
axD.set_yscale("log")
axD.set_xlim(2, 13); axD.set_ylim(0.75, 4.2)
axD.set_xlabel(r"group size $m$")
axD.set_ylabel(r"$\chi\,/\,\chi^{\rm exact}$")
axD.set_xticks([3, 5, 8, 12])
axD.set_yticks([1, 2, 3, 4]); axD.set_yticklabels(["1", "2", "3", "4"])
axD.tick_params(which="both", top=False, right=False)
axD.legend(loc="upper left", frameon=False, fontsize=7.0, handlelength=1.5,
           handletextpad=0.5, labelspacing=0.3, borderaxespad=0.5)
axD.text(0.965, 0.94, "points: Gillespie on the\nmultiplex hypergraph\n(agree with exact, $\\leq$1.3$\\sigma$)",
         transform=axD.transAxes, fontsize=7.4, color=SEC, ha="right",
         va="top", linespacing=1.4)
tag(axD, "d")

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
