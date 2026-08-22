"""Figure 6: two structural questions that only a multiplex can pose.

(a) Channel allocation. Fix the total transmission budget sum_a w_a and move it
    between the layers. The threshold is minimised -- the configuration is most
    dangerous -- at an INTERIOR split in every case tested, symmetric or not:
    spending the whole budget on the single strongest channel is never the worst
    thing one can do, because a pure allocation switches off the cross-layer
    entries of N. Curves are normalised by the better of the two pure
    allocations, so the dip below 1 is exactly what mixing buys the epidemic.

(b, c) Group granularity. "Few large groups" versus "many small groups" has no
    answer until one says what is held fixed, and the two natural choices point
    OPPOSITE ways:
      (b) fixed number of groups per node (k): larger groups are far more
          dangerous, as the superlinear growth of C in m suggests;
      (c) fixed contact budget k(m-1) -- the same number of distinct neighbours
          -- reverses it: lambda_c RISES with m, so many small groups are the
          more dangerous arrangement. The group-level branching factor is
          X*C = (k-1)*C(m,lambda); going to bigger groups gains a superlinear C
          but pays a linear collapse of the excess degree k-1, and the collapse
          wins. At k = 1 the layer cannot spread on its own at all, however
          large its groups.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from scipy.optimize import brentq

from theory import cascade_C, next_gen_matrix, spectral_radius

INK, SEC, MUTED = "#20201e", "#565550", "#9b9a93"
TEAL, CORAL, BLUE, PLUM = "#1C9B8E", "#E76F51", "#2F5FD0", "#7D6B9E"

CASES = [("$m=(3,3)$, $k=(3,3)$", {(3, 3): 1.0}, (3, 3), TEAL),
         ("$m=(2,5)$, $k=(4,2)$", {(4, 2): 1.0}, (2, 5), CORAL),
         ("$m=(3,5)$, $k=(3,3)$", {(3, 3): 1.0}, (3, 5), BLUE),
         ("$m=(3,3)$, $k=(4,2)$", {(4, 2): 1.0}, (3, 3), PLUM)]

BUDGET = 12          # contact budget k(m-1) for panel (c)
KFIX = 4             # groups per node for panel (b)


def lc_joint(P, m, w):
    return brentq(lambda L: spectral_radius(next_gen_matrix(P, m, L, w=np.array(w))) - 1,
                  1e-9, 400, xtol=1e-14)


def lc_single(m, k):
    X = k - 1.0
    if X <= 0:
        return None
    return brentq(lambda L: cascade_C(m, L, 1) * X - 1.0, 1e-9, 400, xtol=1e-14)


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
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(11.4, 3.5))
    fig.subplots_adjust(left=0.058, right=0.987, bottom=0.155, top=0.9, wspace=0.315)

    def tag(ax, s):
        ax.text(-0.20, 1.045, s, transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="bottom", ha="left", color=INK)

    def minor(ax, nx=2, ny=2):
        ax.xaxis.set_minor_locator(AutoMinorLocator(nx))
        ax.yaxis.set_minor_locator(AutoMinorLocator(ny))
        ax.tick_params(which="both", top=False, right=False)

    # ===================================================== (a) allocation
    ss = np.linspace(0.004, 0.996, 320)
    print("=== (a) channel allocation, budget w1+w2 = 1 ===")
    for lab, P, m, col in CASES:
        v = np.array([lc_joint(P, m, (s, 1 - s)) for s in ss])
        pure = min(lc_joint(P, m, (1e-9, 1)), lc_joint(P, m, (1, 1e-9)))
        y = v / pure
        i = y.argmin()
        axA.plot(ss, y, "-", color=col, lw=1.7, zorder=3, label=lab)
        axA.plot([ss[i]], [y[i]], "o", ms=5.0, mfc=col, mec="white", mew=0.8, zorder=5)
        print(f"  {lab}: worst split w1={ss[i]:.3f}, lambda_c/pure={y[i]:.4f} "
              f"({(1-y[i])*100:.2f}% lower than the best pure allocation)")
    axA.axhline(1.0, color=MUTED, lw=0.9, ls=(0, (4, 2.6)), zorder=2)
    axA.set_xlim(0, 1)
    axA.set_yscale("log")
    axA.set_ylim(0.66, 3.75)          # headroom so the legend clears every curve
    axA.set_yticks([0.7, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0])
    axA.set_yticklabels(["0.7", "0.8", "1.0", "1.5", "2.0", "2.5", "3.0"])
    axA.set_xlabel(r"$w_1/(w_1{+}w_2)$")
    axA.set_ylabel(r"$\lambda_c\,/\,\lambda_c^{\mathrm{pure}}$")
    axA.xaxis.set_minor_locator(AutoMinorLocator(2))
    axA.tick_params(which="both", top=False, right=False)
    axA.legend(loc="upper center", ncol=2, fontsize=7.4, handlelength=1.4,
               handletextpad=0.45, labelspacing=0.35, columnspacing=1.0,
               borderaxespad=0.5)
    tag(axA, "a")

    # ============================== (b) granularity at fixed groups-per-node
    ms = [2, 3, 4, 5, 6, 8]
    lb = [lc_single(m, KFIX) for m in ms]
    print(f"\n=== (b) fixed k={KFIX}: bigger groups more dangerous ===")
    for m, v in zip(ms, lb):
        print(f"  m={m}: lambda_c={v:.6f}")
    axB.plot(ms, lb, "-o", color=CORAL, lw=1.7, ms=5.2, mfc=CORAL,
             mec="white", mew=0.8, zorder=3)
    axB.set_yscale("log")
    axB.set_xlabel(r"$m$")
    axB.set_ylabel(r"$\lambda_c$")
    axB.set_xticks(ms)
    axB.tick_params(which="both", top=False, right=False)
    axB.set_title(rf"$k={KFIX}$", fontsize=10, pad=6)
    tag(axB, "b")

    # ============================== (c) granularity at fixed contact budget
    cfg = [(2, 12), (3, 6), (4, 4), (5, 3), (7, 2)]
    mc = [m for m, _ in cfg]
    lc = [lc_single(m, k) for m, k in cfg]
    print(f"\n=== (c) fixed contact budget k(m-1)={BUDGET}: reversed ===")
    for (m, k), v in zip(cfg, lc):
        print(f"  m={m:>2d}, k={k:>2d}: lambda_c={v:.6f}")
    print(f"  m=13, k= 1: X=0, the layer alone can never spread")
    axC.plot(mc, lc, "-o", color=TEAL, lw=1.7, ms=5.2, mfc=TEAL,
             mec="white", mew=0.8, zorder=3)
    for (m, k), v in zip(cfg, lc):
        axC.annotate(rf"$k={k}$", xy=(m, v), xytext=(7, -4),
                     textcoords="offset points", fontsize=7.8, color=SEC,
                     ha="left", va="top")
    axC.set_xlabel(r"$m$")
    axC.set_ylabel(r"$\lambda_c$")
    axC.set_xticks(mc)
    axC.set_ylim(0.088, 0.162)
    minor(axC, nx=2, ny=2)
    axC.set_title(rf"$k(m{{-}}1)={BUDGET}$", fontsize=10, pad=6)
    tag(axC, "c")

    fig.savefig("figure6_design.pdf", bbox_inches="tight")
    fig.savefig("figure6_design.png", bbox_inches="tight", dpi=210)
    print("\nsaved figure6_design.png/.pdf")


if __name__ == "__main__":
    main()
