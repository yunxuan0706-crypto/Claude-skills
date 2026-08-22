"""Figure 5 (core result): the synergy region in the (lambda_1, lambda_2) plane.

With a per-layer rate the threshold condition rho(N) = 1 is no longer a point
but a curve in the plane of the two channel rates. Setting it against the two
single-layer conditions N_aa = 1 -- the rate at which each layer would ignite on
its own -- carves out a wedge:

    lambda_1 < lambda_c^(1),  lambda_2 < lambda_c^(2),  but  rho(N) > 1,

i.e. both channels are individually subcritical yet their union is supercritical.
That wedge is the object of interest: it has no single-layer counterpart, and it
is where a multiplex outbreak cannot be attributed to any one channel.

(a) the anchor family P = {(2,2), (3,3)} at m = (3,3), whose corner point
    lambda_1 = lambda_2 = 0.13 reproduces the N_11 = 0.3860, rho = 1.0132 quoted
    in the text;
(b) the same construction across the inter-layer participation correlation
    rho_12 of Figure 4: the marginals -- hence both single-layer lines -- are
    pinned, so only the critical curve moves, and the wedge widens with rho_12.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from scipy.optimize import brentq

from theory import next_gen_matrix, spectral_radius, cascade_C, degree_moments

# ---------------------------------------------------------------- palette
INK, SEC, MUTED = "#20201e", "#565550", "#9b9a93"
TEAL, CORAL, BLUE = "#1C9B8E", "#E76F51", "#2F5FD0"
FILL_SYN = "#f6d9cc"          # synergy wedge (coral family)
FILL_SUB = "#eef2f6"          # jointly subcritical (cool)

P_ANCHOR = {(2, 2): 0.5, (3, 3): 0.5}
M_GROUPS = (3, 3)


def rho2(P, m, l1, l2):
    return spectral_radius(next_gen_matrix(P, m, 1.0, w=np.array([l1, l2])))


def lc_alone(P, m, a):
    """Rate at which layer a would be critical on its own: N_aa = 1."""
    mk, cross = degree_moments(P, len(m))
    X = cross[a, a] / mk[a] - 1.0
    return brentq(lambda L: cascade_C(m[a], L, 1) * X - 1.0, 1e-9, 60, xtol=1e-14)


def critical_curve(P, m, l1s):
    """lambda_2*(lambda_1) solving rho(N) = 1."""
    out = []
    for l1 in l1s:
        try:
            out.append(brentq(lambda x: rho2(P, m, l1, x) - 1.0,
                              1e-12, 80.0, xtol=1e-14))
        except ValueError:
            out.append(np.nan)
    return np.array(out)


def Pfam(rho):
    """Fixed {2,4} marginals, Pearson correlation rho (the Figure-4 family)."""
    tab = [((2, 2), (1 + rho) / 4), ((4, 4), (1 + rho) / 4),
           ((2, 4), (1 - rho) / 4), ((4, 2), (1 - rho) / 4)]
    return {k: p for k, p in tab if p > 1e-12}


def main():
    mpl.rcParams.update({
        "figure.dpi": 140, "savefig.dpi": 340,
        "font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 9.5,
        "axes.edgecolor": "#c6c5bf", "axes.linewidth": 0.7, "axes.labelcolor": INK,
        "axes.labelpad": 3.5, "text.color": INK,
        "xtick.color": SEC, "ytick.color": SEC,
        "xtick.major.size": 3.4, "ytick.major.size": 3.4,
        "xtick.minor.size": 1.9, "ytick.minor.size": 1.9,
        "axes.grid": True, "grid.color": "#ecebe5", "grid.linewidth": 0.55,
        "axes.axisbelow": True, "legend.frameon": False,
        "lines.solid_capstyle": "round",
    })
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.7, 4.0))
    fig.subplots_adjust(left=0.072, right=0.985, bottom=0.125, top=0.93, wspace=0.26)

    def tag(ax, s):
        ax.text(-0.155, 1.035, s, transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="bottom", ha="left", color=INK)

    def minorticks(ax):
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.tick_params(which="both", top=False, right=False)

    # ================================================== (a) the phase diagram
    lc1 = lc_alone(P_ANCHOR, M_GROUPS, 0)
    lc2 = lc_alone(P_ANCHOR, M_GROUPS, 1)
    hi = 0.46
    l1 = np.linspace(1e-4, lc1 * 0.9995, 400)
    l2c = critical_curve(P_ANCHOR, M_GROUPS, l1)

    # shade: jointly subcritical (below the curve) and the synergy wedge
    axA.fill_between(l1, 0, l2c, color=FILL_SUB, lw=0, zorder=0)
    axA.fill_between(l1, l2c, lc2, where=l2c < lc2, color=FILL_SYN, lw=0, zorder=1)

    axA.axvline(lc1, color=MUTED, lw=1.0, ls=(0, (5, 3)), zorder=2)
    axA.axhline(lc2, color=MUTED, lw=1.0, ls=(0, (5, 3)), zorder=2)
    axA.plot(l1, l2c, "-", color=BLUE, lw=2.0, zorder=4,
             label=r"$\rho(N)=1$")

    lcj = brentq(lambda L: rho2(P_ANCHOR, M_GROUPS, L, L) - 1.0, 1e-9, 9, xtol=1e-14)
    axA.plot([0, hi], [0, hi], "-", color=MUTED, lw=0.8, alpha=0.8, zorder=3)
    axA.plot([lcj], [lcj], "o", ms=6.0, mfc="white", mec=BLUE, mew=1.4, zorder=6)
    axA.plot([0.13], [0.13], "D", ms=5.6, mfc=CORAL, mec="white", mew=0.8, zorder=6)
    axA.annotate(r"$\lambda_1=\lambda_2=0.13$" "\n" r"$\rho(N)=1.0132$",
                 xy=(0.13, 0.13), xytext=(0.207, 0.132), fontsize=8.2, color=SEC,
                 ha="left", va="center", linespacing=1.35,
                 arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7,
                                 shrinkA=0, shrinkB=3))
    axA.text(0.255, 0.335, r"$\rho(N)>1$", fontsize=11, color="#8c3a18",
             ha="center", va="center")
    axA.text(0.052, 0.038, r"$\rho(N)<1$", fontsize=10, color="#46708f",
             ha="left", va="bottom")
    axA.text(lc1 + 0.008, 0.017, r"$\lambda_c^{(1)}=%.4f$" % lc1, fontsize=8.0,
             color=SEC, rotation=90, ha="left", va="bottom")
    axA.text(0.017, lc2 + 0.008, r"$\lambda_c^{(2)}=%.4f$" % lc2, fontsize=8.0,
             color=SEC, ha="left", va="bottom")
    axA.set_xlim(0, hi); axA.set_ylim(0, hi); axA.set_aspect("equal")
    axA.set_xlabel(r"layer-1 rate  $\lambda_1$")
    axA.set_ylabel(r"layer-2 rate  $\lambda_2$")
    minorticks(axA)
    axA.legend(loc="upper right", fontsize=8.8, handlelength=1.7,
               handletextpad=0.6, borderaxespad=0.8)
    tag(axA, "a")

    # ============================== (b) the wedge across the Figure-4 family
    lcA = lc_alone(Pfam(0.0), M_GROUPS, 0)          # marginals pinned -> same for all rho
    hi2 = 0.30
    l1b = np.linspace(1e-4, lcA * 0.9995, 400)
    shades = [(-1.0, TEAL, r"$\rho_{12}=-1$"),
              (0.0, BLUE, r"$\rho_{12}=0$"),
              (1.0, CORAL, r"$\rho_{12}=+1$")]
    curves = {}
    for r, col, lab in shades:
        c = critical_curve(Pfam(r), M_GROUPS, l1b)
        curves[r] = c
        axB.plot(l1b, c, "-", color=col, lw=1.9, zorder=4, label=lab)
    # the widening of the wedge: between rho=-1 and rho=+1
    axB.fill_between(l1b, curves[1.0], curves[-1.0], color="#f2efe6", lw=0, zorder=1)
    axB.axvline(lcA, color=MUTED, lw=1.0, ls=(0, (5, 3)), zorder=2)
    axB.axhline(lcA, color=MUTED, lw=1.0, ls=(0, (5, 3)), zorder=2)
    axB.text(lcA + 0.006, 0.012, r"$\lambda_c^{(a)}$ alone $=%.4f$" % lcA,
             fontsize=8.0, color=SEC, rotation=90, ha="left", va="bottom")
    axB.set_xlim(0, hi2); axB.set_ylim(0, hi2); axB.set_aspect("equal")
    axB.set_xlabel(r"layer-1 rate  $\lambda_1$")
    axB.set_ylabel(r"layer-2 rate  $\lambda_2$")
    minorticks(axB)
    axB.legend(loc="lower left", fontsize=8.6, handlelength=1.7,
               handletextpad=0.6, labelspacing=0.45, borderaxespad=0.9)
    tag(axB, "b")

    fig.savefig("figure5_synergy.pdf", bbox_inches="tight")
    fig.savefig("figure5_synergy.png", bbox_inches="tight", dpi=210)

    # ---- numbers quoted in the text ----
    print(f"anchor family: lc_alone = {lc1:.6f}, joint (l1=l2) = {lcj:.6f}, "
          f"ratio = {lcj/lc1:.4f}  ({(1-lcj/lc1)*100:.1f}% below)")
    print(f"rho(N) at (0.13,0.13) = {rho2(P_ANCHOR,M_GROUPS,0.13,0.13):.6f}")
    for r, _c, _l in shades:
        lj = brentq(lambda L: rho2(Pfam(r), M_GROUPS, L, L) - 1.0, 1e-9, 9, xtol=1e-14)
        print(f"  rho12={r:+.0f}: alone={lcA:.6f}  joint={lj:.6f}  ratio={lj/lcA:.4f}")
    print("saved figure5_synergy.png/.pdf")


if __name__ == "__main__":
    main()
