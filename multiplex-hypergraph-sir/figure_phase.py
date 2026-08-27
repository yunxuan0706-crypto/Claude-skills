"""Showcase: the threshold as a phase portrait -- heatmap + white contours in
the style of the reference boards. Pure theory (spectral radius of the
next-generation matrix), so every pixel is exact and cheap.

Consolidates the physics of the old figures 4 and 5 into one figure:
  (a) rho(N) over the (lambda_1, lambda_2) plane for the anchor family, with the
      rho(N)=1 critical curve, the two single-layer lines N_aa=1, and the
      synergy wedge between them;
  (b) the same surface for the {2,4} correlation family at rho12 = +1, showing
      how correlation slides the critical curve inward;
  (c) the phase boundary rho(N)=1 for rho12 = -1, 0, +1 overlaid;
  (d) Delta lambda_2,c(lambda_1, rho12): how far correlation moves the critical
      layer-2 rate, relative to the uncorrelated case, as a 2D surface.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from scipy.optimize import brentq

from theory import next_gen_matrix, spectral_radius, cascade_C, degree_moments

INK, SEC, MUTED = "#20201e", "#565550", "#9b9a93"
TEAL, CORAL, BLUE = "#1C9B8E", "#E76F51", "#2F5FD0"

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


def rho_grid(P, ax_vals):
    Z = np.empty((len(ax_vals), len(ax_vals)))
    for iy, l2 in enumerate(ax_vals):
        for ix, l1 in enumerate(ax_vals):
            Z[iy, ix] = rho2(P, l1, l2)
    return Z


def main():
    mpl.rcParams.update({
        "figure.dpi": 140, "savefig.dpi": 340,
        "font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 10,
        "axes.edgecolor": "#c6c5bf", "axes.linewidth": 0.7, "axes.labelcolor": INK,
        "axes.labelpad": 3.5, "text.color": INK,
        "xtick.color": SEC, "ytick.color": SEC,
        "xtick.major.size": 3.2, "ytick.major.size": 3.2,
        "axes.axisbelow": True, "legend.frameon": False,
        "lines.solid_capstyle": "round",
    })
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 8.0))
    (axA, axB), (axC, axD) = axes
    fig.subplots_adjust(left=0.085, right=0.965, bottom=0.075, top=0.945,
                        wspace=0.30, hspace=0.28)

    def tag(ax, s, x=-0.16):
        ax.text(x, 1.04, s, transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="bottom", ha="left", color=INK)

    # ------------------------------------------------ (a) anchor phase portrait
    hi = 0.46
    g = np.linspace(0.001, hi, 260)
    Za = rho_grid(P_ANCHOR, g)
    lc1 = lc_alone(P_ANCHOR, 0)
    lc2 = lc_alone(P_ANCHOR, 1)
    l1 = np.linspace(1e-4, lc1 * 0.9995, 320)
    curve = crit_curve(P_ANCHOR, l1)

    im = axA.pcolormesh(g, g, Za, cmap="magma", vmin=0, vmax=3.0,
                        shading="gouraud", rasterized=True)
    cs = axA.contour(g, g, Za, levels=[0.25, 0.5, 1.0, 1.5, 2.0, 2.5],
                     colors="white", linewidths=0.7, alpha=0.8)
    axA.clabel(cs, [0.5, 1.5, 2.5], fmt="%.1f", fontsize=7, colors="white")
    axA.plot(l1, curve, color="white", lw=2.4)
    axA.plot(l1, curve, color=INK, lw=1.0, ls=(0, (1, 1.4)))
    axA.axvline(lc1, color="white", lw=0.9, ls=(0, (4, 3)), alpha=0.7)
    axA.axhline(lc2, color="white", lw=0.9, ls=(0, (4, 3)), alpha=0.7)
    axA.plot([0.13], [0.13], "o", ms=6, mfc="none", mec="white", mew=1.6)
    axA.text(0.135, 0.115, r"$\rho{=}1.0132$", color="white", fontsize=8,
             ha="left", va="top")
    axA.text(0.30, 0.33, "synergy", color="white", fontsize=10.5, style="italic",
             ha="center", rotation=-38)
    axA.set_xlim(0, hi); axA.set_ylim(0, hi); axA.set_aspect("equal")
    axA.set_xlabel(r"$\lambda_1$"); axA.set_ylabel(r"$\lambda_2$")
    axA.set_title(r"$P=\{(2,2),(3,3)\}$", fontsize=10, pad=5)
    cb = fig.colorbar(im, ax=axA, fraction=0.046, pad=0.03)
    cb.set_label(r"$\rho(N)$", fontsize=9)
    cb.ax.tick_params(labelsize=7.5)
    tag(axA, "a")

    # ------------------------------- (b) correlated family, rho12 = +1
    hi2 = 0.30
    g2 = np.linspace(0.001, hi2, 260)
    Zb = rho_grid(Pf(1.0), g2)
    lcA = lc_alone(Pf(0.0), 0)
    l1b = np.linspace(1e-4, lcA * 0.9995, 320)
    curveb = crit_curve(Pf(1.0), l1b)
    im2 = axB.pcolormesh(g2, g2, Zb, cmap="magma", vmin=0, vmax=3.0,
                         shading="gouraud", rasterized=True)
    cs2 = axB.contour(g2, g2, Zb, levels=[0.25, 0.5, 1.0, 1.5, 2.0, 2.5],
                      colors="white", linewidths=0.7, alpha=0.8)
    axB.clabel(cs2, [0.5, 1.5, 2.5], fmt="%.1f", fontsize=7, colors="white")
    axB.plot(l1b, curveb, color="white", lw=2.4)
    axB.plot(l1b, curveb, color=INK, lw=1.0, ls=(0, (1, 1.4)))
    axB.axvline(lcA, color="white", lw=0.9, ls=(0, (4, 3)), alpha=0.7)
    axB.axhline(lcA, color="white", lw=0.9, ls=(0, (4, 3)), alpha=0.7)
    axB.set_xlim(0, hi2); axB.set_ylim(0, hi2); axB.set_aspect("equal")
    axB.set_xlabel(r"$\lambda_1$"); axB.set_ylabel(r"$\lambda_2$")
    axB.set_title(r"$\{2,4\}$ family,  $\rho_{12}=+1$", fontsize=10, pad=5)
    cb2 = fig.colorbar(im2, ax=axB, fraction=0.046, pad=0.03)
    cb2.set_label(r"$\rho(N)$", fontsize=9)
    cb2.ax.tick_params(labelsize=7.5)
    tag(axB, "b")

    # ------------------------------- (c) phase boundary shift by rho12
    for r, col, lab in [(-1.0, TEAL, r"$\rho_{12}=-1$"),
                        (0.0, BLUE, r"$\rho_{12}=0$"),
                        (1.0, CORAL, r"$\rho_{12}=+1$")]:
        axC.plot(l1b, crit_curve(Pf(r), l1b), color=col, lw=2.0, label=lab)
    axC.axvline(lcA, color=MUTED, lw=0.9, ls=(0, (4, 3)))
    axC.axhline(lcA, color=MUTED, lw=0.9, ls=(0, (4, 3)))
    axC.set_xlim(0, hi2); axC.set_ylim(0, hi2); axC.set_aspect("equal")
    axC.set_xlabel(r"$\lambda_1$"); axC.set_ylabel(r"$\lambda_2$")
    axC.set_title(r"phase boundary $\rho(N)=1$", fontsize=10, pad=5)
    axC.xaxis.set_minor_locator(AutoMinorLocator(2))
    axC.yaxis.set_minor_locator(AutoMinorLocator(2))
    axC.tick_params(which="both", top=False, right=False)
    axC.grid(True, color="#ecebe5", lw=0.55)
    axC.legend(loc="upper right", fontsize=8.6, handlelength=1.6,
               handletextpad=0.6, labelspacing=0.45)
    tag(axC, "c")

    # ------------------------------- (d) Delta lambda2,c(lambda1, rho12) surface
    rr = np.linspace(-1, 1, 120)
    l1d = np.linspace(0.002, lcA * 0.985, 120)
    base = crit_curve(Pf(0.0), l1d)
    D = np.empty((len(rr), len(l1d)))
    for iy, r in enumerate(rr):
        cr = crit_curve(Pf(r), l1d)
        D[iy, :] = cr - base
    vmax = np.nanmax(np.abs(D))
    imd = axD.pcolormesh(l1d, rr, D, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                         shading="gouraud", rasterized=True)
    csd = axD.contour(l1d, rr, D, levels=[-0.02, -0.01, 0.0, 0.01, 0.02],
                      colors="k", linewidths=0.5, alpha=0.5)
    axD.clabel(csd, fmt="%.2f", fontsize=6.5)
    axD.set_xlim(0, lcA * 0.985); axD.set_ylim(-1, 1)
    axD.set_xlabel(r"$\lambda_1$"); axD.set_ylabel(r"$\rho_{12}$")
    axD.set_title(r"$\Delta\lambda_{2,c}$  (rel. to $\rho_{12}{=}0$)", fontsize=10, pad=5)
    cbd = fig.colorbar(imd, ax=axD, fraction=0.046, pad=0.03)
    cbd.set_label(r"$\Delta\lambda_{2,c}$", fontsize=9)
    cbd.ax.tick_params(labelsize=7.5)
    tag(axD, "d")

    fig.savefig("figure_phase.pdf", bbox_inches="tight")
    fig.savefig("figure_phase.png", bbox_inches="tight", dpi=200)
    print("saved figure_phase.png/.pdf")
    print(f"anchor: lc_alone={lc1:.5f}, joint={crit_curve(P_ANCHOR,[lc1*0.99999])[0]:.5f}")
    print(f"{'{2,4}'}: lc_alone={lcA:.5f}; Delta range [{D.min():+.4f},{D.max():+.4f}]")


if __name__ == "__main__":
    main()
