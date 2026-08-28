"""Figure 5: structural design, and where the theory stops.

  (a) channel allocation. At fixed total budget the most dangerous split is
      interior in every case tested -- a pure allocation zeroes the cross-layer
      entries of N. Curves normalised by the better pure allocation.
  (b) group granularity. "Few large groups vs many small groups" reverses under
      the two natural normalisations: fixed groups-per-node k makes big groups
      more dangerous; fixed contact budget k(m-1) makes them less.
  (c) inter-layer overlap o -- which N cannot see -- moves lambda_c anyway:
      falsifiability of the tree closure. Measured lambda_c (scaled errors)
      against the overlap-blind prediction and the o=1 analytic limit.
  (d) mechanism: the group branching factor the theory counts, 3C(3,lambda),
      against what the o=1 process delivers, C(3,2lambda); each crosses unity at
      its own threshold.

Design panels (a,b) are pure theory; the overlap data (c) is figure7_data.json.
"""
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from scipy.optimize import brentq

from theory import cascade_C, next_gen_matrix, spectral_radius

INK, SEC, MUTED = "#20201e", "#565550", "#9b9a93"
TEAL, CORAL, BLUE, PLUM = "#1C9B8E", "#E76F51", "#2F5FD0", "#7D6B9E"

CASES = [(r"$m{=}(3,3),\,k{=}(3,3)$", {(3, 3): 1.0}, (3, 3), TEAL),
         (r"$m{=}(2,5),\,k{=}(4,2)$", {(4, 2): 1.0}, (2, 5), CORAL),
         (r"$m{=}(3,5),\,k{=}(3,3)$", {(3, 3): 1.0}, (3, 5), BLUE)]
BUDGET = 12


def lc_joint(P, m, w):
    return brentq(lambda L: spectral_radius(next_gen_matrix(P, m, L, w=np.array(w))) - 1,
                  1e-9, 400, xtol=1e-14)


def lc_single(m, k):
    X = k - 1.0
    return None if X <= 0 else brentq(lambda L: cascade_C(m, L, 1) * X - 1.0, 1e-9, 400, xtol=1e-14)


def wls(x, y, ye):
    W = 1 / ye ** 2; X = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(X.T @ (W[:, None] * X)); a, b = cov @ (X.T @ (W * y)); return a, b


def main():
    mpl.rcParams.update({
        "figure.dpi": 140, "savefig.dpi": 340,
        "font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 9.5,
        "axes.edgecolor": "#c6c5bf", "axes.linewidth": 0.7, "axes.labelcolor": INK,
        "axes.labelpad": 3.0, "text.color": INK,
        "xtick.color": SEC, "ytick.color": SEC,
        "xtick.major.size": 3.2, "ytick.major.size": 3.2,
        "xtick.minor.size": 1.8, "ytick.minor.size": 1.8,
        "axes.grid": True, "grid.color": "#ecebe5", "grid.linewidth": 0.55,
        "axes.axisbelow": True, "legend.frameon": False,
        "lines.solid_capstyle": "round",
    })
    fig, ((axA, axB), (axC, axD)) = plt.subplots(2, 2, figsize=(9.3, 7.0))
    fig.subplots_adjust(left=0.085, right=0.975, bottom=0.085, top=0.955,
                        wspace=0.28, hspace=0.30)

    def tag(ax, s):
        ax.text(-0.155, 1.04, s, transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="bottom", ha="left", color=INK)

    def minor(ax, nx=2, ny=2):
        ax.xaxis.set_minor_locator(AutoMinorLocator(nx))
        ax.yaxis.set_minor_locator(AutoMinorLocator(ny))
        ax.tick_params(which="both", top=False, right=False)

    # -------------------------------------------------- (a) channel allocation
    ss = np.linspace(0.004, 0.996, 320)
    for lab, P, m, col in CASES:
        v = np.array([lc_joint(P, m, (s, 1 - s)) for s in ss])
        pure = min(lc_joint(P, m, (1e-9, 1)), lc_joint(P, m, (1, 1e-9)))
        y = v / pure; i = y.argmin()
        axA.plot(ss, y, "-", color=col, lw=1.6, zorder=3, label=lab)
        axA.plot([ss[i]], [y[i]], "o", ms=4.8, mfc=col, mec="white", mew=0.8, zorder=5)
    axA.axhline(1.0, color=MUTED, lw=0.9, ls=(0, (4, 2.6)), zorder=2)
    axA.set_xlim(0, 1); axA.set_yscale("log"); axA.set_ylim(0.66, 3.75)
    axA.set_yticks([0.7, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0])
    axA.set_yticklabels(["0.7", "0.8", "1.0", "1.5", "2.0", "2.5", "3.0"])
    axA.set_xlabel(r"$w_1/(w_1{+}w_2)$"); axA.set_ylabel(r"$\lambda_c/\lambda_c^{\mathrm{pure}}$")
    axA.xaxis.set_minor_locator(AutoMinorLocator(2)); axA.tick_params(which="both", top=False, right=False)
    axA.legend(loc="upper left", fontsize=7.6, handlelength=1.5,
               handletextpad=0.5, labelspacing=0.4, borderaxespad=0.7)
    tag(axA, "a")

    # -------------------------------------------------- (b) granularity reversal
    msK = [2, 3, 4, 5, 6, 8]
    lbK = [lc_single(m, 4) for m in msK]
    cfgB = [(2, 12), (3, 6), (4, 4), (5, 3), (7, 2)]
    msB = [m for m, _ in cfgB]; lbB = [lc_single(m, k) for m, k in cfgB]
    axB.plot(msK, lbK, "-o", color=CORAL, lw=1.7, ms=5.0, mfc=CORAL, mec="white",
             mew=0.8, zorder=3, label=r"fixed $k=4$")
    axB.plot(msB, lbB, "-s", color=TEAL, lw=1.7, ms=5.0, mfc=TEAL, mec="white",
             mew=0.8, zorder=3, label=r"fixed $k(m{-}1)=12$")
    axB.set_yscale("log")
    axB.set_xlabel(r"group size $m$"); axB.set_ylabel(r"$\lambda_c$  (single layer)")
    axB.set_xticks(range(2, 9)); axB.tick_params(which="both", top=False, right=False)
    axB.legend(loc="lower left", fontsize=8.0, handlelength=1.6, handletextpad=0.5,
               labelspacing=0.45)
    tag(axB, "b")

    # -------------------------------------------------- (c) overlap falsifiability
    d7 = json.load(open("figure7_data.json"))
    o = np.array(d7["o"]); lc = np.array(d7["lc"]); se = np.array(d7["se"])
    fac = []
    for i in range(len(o)):
        x, y, ye = map(np.array, (d7["lams_all"][i], d7["y_all"][i], d7["ye_all"][i]))
        a, b = wls(x, y, ye); chi2 = (((y - (a + b * x)) / ye) ** 2).sum() / (len(x) - 2)
        fac.append(max(1.0, np.sqrt(chi2)))
    se_s = se * np.array(fac)
    axC.axhline(d7["lc_theory"], color=BLUE, lw=1.6, zorder=2)
    axC.axhline(d7["lc_o1"], color=CORAL, lw=1.6, ls=(0, (4, 2.4)), zorder=2)
    axC.errorbar(o, lc, yerr=se_s, fmt="o", ms=5.0, mfc=PLUM, mec="white", mew=0.8,
                 ecolor=PLUM, elinewidth=1.1, capsize=2.4, capthick=0.9, zorder=4)
    axC.plot(o, lc, "-", color=PLUM, lw=1.2, zorder=3)
    axC.set_xlim(-0.06, 1.06); axC.set_ylim(d7["lc_theory"] - 0.02, d7["lc_o1"] + 0.03)
    axC.set_xlabel(r"inter-layer overlap $o$"); axC.set_ylabel(r"$\lambda_c$")
    minor(axC)
    # inline labels on the two reference lines, placed clear of the curve and edges
    axC.text(0.97, d7["lc_theory"] + 0.006, r"$\rho(N)=1$", color=BLUE, fontsize=8.4,
             ha="right", va="bottom")
    axC.text(0.03, d7["lc_o1"] - 0.006, r"$C(3,2\lambda_c)=1$", color=CORAL, fontsize=8.4,
             ha="left", va="top")
    tag(axC, "c")

    # -------------------------------------------------- (d) mechanism
    lam = np.linspace(0.0, 0.55, 260)
    th = np.array([3 * cascade_C(3, L, 1) for L in lam])
    ac = np.array([cascade_C(3, 2 * L, 1) for L in lam])
    axD.axhline(1.0, color=MUTED, lw=0.9, ls=(0, (4, 2.4)), zorder=2)
    axD.plot(lam, th, "-", color=BLUE, lw=1.8, zorder=3, label=r"$3\,C(3,\lambda)$")
    axD.plot(lam, ac, "-", color=CORAL, lw=1.8, zorder=3, label=r"$C(3,2\lambda)$")
    axD.plot([d7["lc_theory"]], [1.0], "o", ms=5.2, mfc="white", mec=BLUE, mew=1.2, zorder=5)
    axD.plot([d7["lc_o1"]], [1.0], "o", ms=5.2, mfc="white", mec=CORAL, mew=1.2, zorder=5)
    axD.set_xlim(0, 0.55); axD.set_ylim(0, 2.6)
    axD.set_xlabel(r"$\lambda$"); axD.set_ylabel(r"group branching factor")
    minor(axD)
    axD.legend(loc="upper left", fontsize=8.6, handlelength=1.7, handletextpad=0.6,
               labelspacing=0.5)
    tag(axD, "d")

    fig.savefig("figure4_structure.pdf", bbox_inches="tight")
    fig.savefig("figure4_structure.png", bbox_inches="tight", dpi=200)
    print("saved figure4_structure.png/.pdf")


if __name__ == "__main__":
    main()
