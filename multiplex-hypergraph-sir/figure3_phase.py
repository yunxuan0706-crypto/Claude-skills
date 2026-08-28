"""Figure 3 (central result): the epidemic threshold as a phase portrait.

With a per-layer rate the condition rho(N)=1 is a curve in the (lambda_1,lambda_2)
plane, not a point. Four panels:
  (a) rho(N) over the plane for the anchor family, magma heatmap + white
      contours, the rho(N)=1 boundary, the single-layer lines N_aa=1, and the
      synergy wedge between them (each layer subcritical alone, together not);
  (b) lambda_c(rho12) for the fixed-marginal {2,4} family: exact rho(N)=1 curve
      vs a bootstrap of the direct Gillespie measurement (cached);
  (c) the rho(N)=1 boundary at rho12 = -1, 0, +1 -- correlation slides it in;
  (d) Delta lambda_2,c(lambda_1, rho12), how far correlation moves the critical
      layer-2 rate relative to the uncorrelated case, as a signed surface.

All theory is the spectral radius of N (exact, cheap); the simulation in (b) is
read from figure4_data.json.
"""
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.patches import Patch
from scipy.optimize import brentq

from theory import next_gen_matrix, spectral_radius, cascade_C, degree_moments, lambda_c_rho
from figure4_rho12 import bootstrap_lc

INK, SEC, MUTED = "#20201e", "#565550", "#9b9a93"
TEAL, CORAL, BLUE = "#1C9B8E", "#E76F51", "#2F5FD0"
MED = "#a23c1a"
P_ANCHOR = {(2, 2): 0.5, (3, 3): 0.5}
M = (3, 3)


def Pf(r):
    tab = [((2, 2), (1 + r) / 4), ((4, 4), (1 + r) / 4),
           ((2, 4), (1 - r) / 4), ((4, 2), (1 - r) / 4)]
    return {k: p for k, p in tab if p > 1e-12}


def rho2(P, l1, l2):
    return spectral_radius(next_gen_matrix(P, M, 1.0, w=np.array([l1, l2])))


def lc_alone(P, a):
    mk, cr = degree_moments(P, 2)
    X = cr[a, a] / mk[a] - 1.0
    return brentq(lambda L: cascade_C(M[a], L, 1) * X - 1.0, 1e-9, 60, xtol=1e-14)


def crit_curve(P, l1s):
    out = []
    for l1 in l1s:
        try:
            out.append(brentq(lambda x: rho2(P, l1, x) - 1.0, 1e-12, 80, xtol=1e-13))
        except ValueError:
            out.append(np.nan)
    return np.array(out)


def rho_grid(P, gv):
    Z = np.empty((len(gv), len(gv)))
    for iy, l2 in enumerate(gv):
        for ix, l1 in enumerate(gv):
            Z[iy, ix] = rho2(P, l1, l2)
    return Z


def main():
    mpl.rcParams.update({
        "figure.dpi": 140, "savefig.dpi": 340,
        "font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 10,
        "axes.edgecolor": "#c6c5bf", "axes.linewidth": 0.7, "axes.labelcolor": INK,
        "axes.labelpad": 3.0, "text.color": INK,
        "xtick.color": SEC, "ytick.color": SEC,
        "xtick.major.size": 3.2, "ytick.major.size": 3.2,
        "axes.axisbelow": True, "legend.frameon": False,
        "lines.solid_capstyle": "round",
    })
    fig, ((axA, axB), (axC, axD)) = plt.subplots(2, 2, figsize=(9.5, 8.0))
    fig.subplots_adjust(left=0.082, right=0.965, bottom=0.072, top=0.955,
                        wspace=0.40, hspace=0.26)

    def tag(ax, s):
        ax.text(-0.15, 1.035, s, transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="bottom", ha="left", color=INK)

    def minor(ax):
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.tick_params(which="both", top=False, right=False)

    # ================================================== (a) anchor phase portrait
    hi = 0.46
    gv = np.linspace(0.001, hi, 260)
    Za = rho_grid(P_ANCHOR, gv)
    lc1, lc2 = lc_alone(P_ANCHOR, 0), lc_alone(P_ANCHOR, 1)
    l1 = np.linspace(1e-4, lc1 * 0.9995, 320)
    cur = crit_curve(P_ANCHOR, l1)
    im = axA.pcolormesh(gv, gv, Za, cmap="magma", vmin=0, vmax=3.0,
                        shading="gouraud", rasterized=True)
    cs = axA.contour(gv, gv, Za, levels=[0.25, 0.5, 1.0, 1.5, 2.0, 2.5],
                     colors="white", linewidths=0.7, alpha=0.85)
    axA.clabel(cs, [0.5, 1.5, 2.5], fmt="%.1f", fontsize=7, colors="white")
    # the synergy wedge S of (3.1): above the critical curve, yet below BOTH
    # single-layer lines -- each channel subcritical alone, supercritical together.
    # Shaded so the central claim of 3.1 is legible without tracing contours.
    axA.fill_between(l1, cur, lc2, where=(cur < lc2), color="white", alpha=0.16,
                     lw=0, zorder=2)
    axA.fill_between(l1, cur, lc2, where=(cur < lc2), facecolor="none",
                     hatch="///", edgecolor="white", lw=0.0, alpha=0.5, zorder=2)
    axA.plot(l1, cur, color="white", lw=2.4, zorder=3)
    axA.plot(l1, cur, color=INK, lw=1.0, ls=(0, (1, 1.4)), zorder=3)
    axA.axvline(lc1, color="white", lw=0.9, ls=(0, (4, 3)), alpha=0.7, zorder=3)
    axA.axhline(lc2, color="white", lw=0.9, ls=(0, (4, 3)), alpha=0.7, zorder=3)
    axA.plot([0.13], [0.13], "o", ms=6, mfc="none", mec="white", mew=1.6, zorder=4)
    axA.text(0.253, 0.253, r"$\mathcal{S}$", color="white", fontsize=15,
             ha="center", va="center", zorder=4)
    axA.text(0.253, 0.196, "synergy", color="white", fontsize=8.6, style="italic",
             ha="center", va="center", zorder=4)
    axA.set_xlim(0, hi); axA.set_ylim(0, hi); axA.set_aspect("equal")
    axA.set_xlabel(r"$\lambda_1$"); axA.set_ylabel(r"$\lambda_2$")
    axA.set_title(r"$P=\{(2,2),(3,3)\}$", fontsize=9.6, pad=4)
    cb = fig.colorbar(im, ax=axA, fraction=0.046, pad=0.03)
    cb.set_label(r"$\rho(N)$", fontsize=9); cb.ax.tick_params(labelsize=7.5)
    tag(axA, "a")

    # ==================================== (b) lambda_c(rho12): theory vs bootstrap
    d4 = json.load(open("figure4_data.json"))
    rhos = np.linspace(-1, 1, 41)
    lc_th = np.array([lambda_c_rho(Pf(r), M) for r in rhos])
    xs = np.array(d4["rho_sim"])
    lams_all = [np.array(v) for v in d4["lams_all"]]
    y_all = [np.array(v) for v in d4["y_all"]]
    ye_all = [np.array(v) for v in d4["ye_all"]]
    rng = np.random.default_rng(2024)
    boots = [bootstrap_lc(lams_all[i], y_all[i], ye_all[i], B=4000, rng=rng)
             for i in range(len(xs))]
    axB.plot(rhos, lc_th, "-", color=BLUE, lw=1.9, zorder=2, label=r"theory $\rho(N)=1$")
    bp = axB.boxplot(boots, positions=xs, widths=0.11, patch_artist=True,
                     whis=(5, 95), showfliers=False, showcaps=True,
                     manage_ticks=False, zorder=4,
                     medianprops=dict(color=MED, lw=1.3),
                     boxprops=dict(lw=0.8, edgecolor=MED),
                     whiskerprops=dict(color=MED, lw=0.8),
                     capprops=dict(color=MED, lw=0.8))
    for p in bp["boxes"]:
        p.set_facecolor("#f4c9b6"); p.set_alpha(0.9)
    axB.set_xlim(-1.15, 1.15); axB.set_xticks([-1, -0.5, 0, 0.5, 1])
    axB.set_xlabel(r"$\rho_{12}$"); axB.set_ylabel(r"$\lambda_c$")
    minor(axB); axB.grid(True, color="#ecebe5", lw=0.55)
    axB.legend([axB.lines[0], Patch(facecolor="#f4c9b6", edgecolor=MED, lw=0.8)],
               [r"theory $\rho(N)=1$", "simulation"], loc="upper right",
               fontsize=8.4, handlelength=1.6, handletextpad=0.6, labelspacing=0.5)
    tag(axB, "b")

    # ==================================== (c) phase boundary shift by rho12
    lcA = lc_alone(Pf(0.0), 0)
    hi2 = 0.30
    l1b = np.linspace(1e-4, lcA * 0.9995, 320)
    for r, col, lab in [(-1.0, TEAL, r"$\rho_{12}=-1$"),
                        (0.0, BLUE, r"$\rho_{12}=0$"),
                        (1.0, CORAL, r"$\rho_{12}=+1$")]:
        axC.plot(l1b, crit_curve(Pf(r), l1b), color=col, lw=2.0, label=lab)
    axC.axvline(lcA, color=MUTED, lw=0.9, ls=(0, (4, 3)))
    axC.axhline(lcA, color=MUTED, lw=0.9, ls=(0, (4, 3)))
    axC.set_xlim(0, hi2); axC.set_ylim(0, hi2); axC.set_aspect("equal")
    axC.set_xlabel(r"$\lambda_1$"); axC.set_ylabel(r"$\lambda_2$")
    axC.set_title(r"phase boundary $\rho(N)=1$", fontsize=9.6, pad=4)
    minor(axC); axC.grid(True, color="#ecebe5", lw=0.55)
    axC.legend(loc="upper right", fontsize=8.4, handlelength=1.6,
               handletextpad=0.6, labelspacing=0.45)
    tag(axC, "c")

    # ==================================== (d) Delta lambda2,c(lambda1, rho12)
    rr = np.linspace(-1, 1, 120)
    l1d = np.linspace(0.002, lcA * 0.985, 120)
    base = crit_curve(Pf(0.0), l1d)
    D = np.array([crit_curve(Pf(r), l1d) - base for r in rr])
    vmax = np.nanmax(np.abs(D))
    imd = axD.pcolormesh(l1d, rr, D, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                         shading="gouraud", rasterized=True)
    csd = axD.contour(l1d, rr, D, levels=[-0.01, 0.0, 0.01], colors="k",
                      linewidths=0.5, alpha=0.45)
    axD.clabel(csd, fmt="%.2f", fontsize=6.5)
    axD.set_xlim(0, lcA * 0.985); axD.set_ylim(-1, 1)
    axD.set_xlabel(r"$\lambda_1$"); axD.set_ylabel(r"$\rho_{12}$")
    axD.set_title(r"$\Delta\lambda_{2,c}$  (rel. to $\rho_{12}{=}0$)", fontsize=9.6, pad=4)
    cbd = fig.colorbar(imd, ax=axD, fraction=0.046, pad=0.03)
    cbd.set_label(r"$\Delta\lambda_{2,c}$", fontsize=9); cbd.ax.tick_params(labelsize=7.5)
    tag(axD, "d")

    fig.savefig("figure3_phase.pdf", bbox_inches="tight")
    fig.savefig("figure3_phase.png", bbox_inches="tight", dpi=200)
    joint = crit_curve(P_ANCHOR, [lc1 * 0.99999])[0]
    print("saved figure3_phase.png/.pdf")
    print(f"anchor lc_alone={lc1:.5f} joint={lambda_c_rho(P_ANCHOR,M):.5f} "
          f"({(1-lambda_c_rho(P_ANCHOR,M)/lc1)*100:.1f}% below)")
    print(f"{{2,4}} lc_alone={lcA:.5f}; theory drop {(1-lc_th[-1]/lc_th[0])*100:.1f}%")


if __name__ == "__main__":
    main()
