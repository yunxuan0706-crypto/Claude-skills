"""Figure 7: the falsifiable test -- does the threshold move with inter-layer
group overlap, which the next-generation matrix cannot see?

(2.7) reads only P(k), m and theta. Two multiplexes with the same P(k) and the
same m therefore carry the same predicted lambda_c, no matter how their layers
are positioned relative to one another. Overlap -- a layer-2 group re-using a
node pair that already shares a layer-1 group -- closes a 4-cycle, exactly what
the local-tree assumption discards. So a measured dependence of lambda_c on o is
a direct, quantitative reading of the tree closure's failure.

Family: every node has layer-degree (2,2), m = (3,3), and a fraction o of the
layer-1 groups is reproduced verbatim in layer 2 (see overlap.py). P(k) and m
are identical throughout; o is the only thing that varies.

Two reference lines bracket the effect:
  theory   1 * C(3,lam) * 3 = 1                  -> lambda_c = 0.18614  (flat in o)
  o = 1    layer 2 is a copy, so each physical group carries rate 2*lam and a
           node has 2 physical groups: 1 * C(3, 2 lam) = 1 -> lambda_c = 0.40974
"""
import json
import os
import numpy as np
from scipy.optimize import brentq

from theory import cascade_C, lambda_c_rho
from overlap import build_overlap_multiplex
from simulate import _membership, outbreak_size
from figure4_rho12 import _wls, _xintercept, bootstrap_lc

M_GROUPS = (3, 3)
KDEG = 2
N_SIM = 20000
NGRAPHS = 4
# The window must stay clear of lambda_c: 1/chi is linear only asymptotically,
# and close to the threshold chi grows until the outbreak cap truncates it. A
# first attempt with [0.72,0.92] put the top point within 3% of lambda_c for the
# intermediate overlaps -- where the initial guess was weakest -- and the fits
# came back with 15-sigma residual curvature. Backing the window off and letting
# it re-centre three times keeps every point comfortably subcritical.
WINDOW = np.array([0.60, 0.66, 0.72, 0.78, 0.84])
O_SIM = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

P_FAM = {(KDEG, KDEG): 1.0}
LC_THEORY = lambda_c_rho(P_FAM, M_GROUPS)                       # 0.186141
LC_O1 = brentq(lambda L: cascade_C(3, L, 1) - 1.0, 1e-9, 60, xtol=1e-14) / 2


def chi_at(o, lam, n_seeds, seed, cap=None):
    """Mean subcritical outbreak size, averaged over NGRAPHS realisations."""
    rng = np.random.default_rng(seed)
    cap = cap or N_SIM // 4          # generous: truncating the tail biases chi low
    sizes = []
    for _ in range(NGRAPHS):
        groups, _o = build_overlap_multiplex(N_SIM, M_GROUPS, KDEG, o, rng)
        owner = _membership(groups, N_SIM)
        beta = [lam, lam]
        for _ in range(n_seeds // NGRAPHS):
            v = int(rng.integers(N_SIM))
            s, _h = outbreak_size(groups, owner, beta, [1, 1], v, rng, cap=cap)
            sizes.append(s)
    a = np.asarray(sizes, float)
    return a.mean(), a.std(ddof=1) / np.sqrt(len(a))


def realised_overlap(o, seed, reps=3):
    rng = np.random.default_rng(seed)
    vals = [build_overlap_multiplex(N_SIM, M_GROUPS, KDEG, o, rng)[1]
            for _ in range(reps)]
    return float(np.mean(vals))


def lambda_c_of_o(o, seed):
    """Three passes: two cheap ones to re-centre the window on lambda_c, then a
    high-statistics measurement. Returns the fit and its residual pulls, so a
    badly placed window cannot pass silently."""
    guess = LC_THEORY + o * (LC_O1 - LC_THEORY)          # bracketing interpolation
    for npass, nseed in ((0, 1200), (1, 2000), (2, 6000)):
        lams = WINDOW * guess
        y, ye = [], []
        for L in lams:
            c, se = chi_at(o, L, nseed, seed + 17 * npass + int(1000 * L))
            y.append(1.0 / c)
            ye.append(se / c ** 2)
        y, ye = np.array(y), np.array(ye)
        lc, selc = _xintercept(lams, y, ye)
        guess = lc
    a, b, _ = _wls(lams, y, ye)
    pulls = (y - (a + b * lams)) / ye
    return lc, selc, lams, y, ye, pulls


def main():
    cache = "figure7_data.json"
    if os.path.exists(cache):
        d = json.load(open(cache))
        os_, lcs, ses = (np.array(d["o"]), np.array(d["lc"]), np.array(d["se"]))
        opair = np.array(d["o_pair"])
        lams_all = [np.array(v) for v in d["lams_all"]]
        y_all = [np.array(v) for v in d["y_all"]]
        ye_all = [np.array(v) for v in d["ye_all"]]
    else:
        os_, lcs, ses, opair = [], [], [], []
        lams_all, y_all, ye_all = [], [], []
        for i, o in enumerate(O_SIM):
            lc, se, lams, y, ye, pulls = lambda_c_of_o(o, seed=500 + 31 * i)
            op = realised_overlap(o, seed=9000 + i)
            os_.append(o); lcs.append(lc); ses.append(se); opair.append(op)
            lams_all.append(lams); y_all.append(y); ye_all.append(ye)
            print(f"o={o:.2f} (pair {op:.4f}): lambda_c={lc:.5f}+-{se:.5f}"
                  f"  theory={LC_THEORY:.5f}  ratio={lc/LC_THEORY:.3f}"
                  f"  | fit max|pull|={np.abs(pulls).max():.1f}"
                  f"  top lam/lc={lams[-1]/lc:.3f}", flush=True)
        os_, lcs, ses, opair = map(np.array, (os_, lcs, ses, opair))
        json.dump({"o": os_.tolist(), "o_pair": opair.tolist(),
                   "lc": lcs.tolist(), "se": ses.tolist(),
                   "lams_all": [v.tolist() for v in lams_all],
                   "y_all": [v.tolist() for v in y_all],
                   "ye_all": [v.tolist() for v in ye_all],
                   "N": N_SIM, "k": KDEG, "m": list(M_GROUPS),
                   "lc_theory": LC_THEORY, "lc_o1": LC_O1},
                  open(cache, "w"), indent=2)

    print(f"\nlambda_c: {lcs[0]:.5f} (o=0) -> {lcs[-1]:.5f} (o=1), "
          f"a {(lcs[-1]/lcs[0]-1)*100:.1f}% rise")
    print(f"theory (flat) = {LC_THEORY:.5f};  o=1 analytic = {LC_O1:.5f}")
    print(f"pull at o=0 : {(lcs[0]-LC_THEORY)/ses[0]:+.2f} sigma  "
          f"(control: the tree closure should be exact here)")
    print(f"pull at o=1 : {(lcs[-1]-LC_THEORY)/ses[-1]:+.2f} sigma  vs the flat prediction")
    plot(os_, opair, lcs, ses, lams_all, y_all, ye_all)
    return os_, opair, lcs, ses, lams_all, y_all, ye_all


# --------------------------------------------------------------------------
def plot(os_, opair, lcs, ses, lams_all, y_all, ye_all):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.cm import ScalarMappable
    from matplotlib.ticker import AutoMinorLocator

    INK, SEC, MUTED = "#20201e", "#565550", "#9b9a93"
    TEAL, CORAL, BLUE = "#1C9B8E", "#E76F51", "#2F5FD0"
    SEQ = LinearSegmentedColormap.from_list("o", ["#dfeae4", TEAL, "#0d4f47"])
    NORM = Normalize(0.0, 1.0)

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
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(11.4, 3.55))
    fig.subplots_adjust(left=0.058, right=0.987, bottom=0.155, top=0.9, wspace=0.315)

    def tag(ax, s):
        ax.text(-0.20, 1.045, s, transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="bottom", ha="left", color=INK)

    def minor(ax, nx=2, ny=2):
        ax.xaxis.set_minor_locator(AutoMinorLocator(nx))
        ax.yaxis.set_minor_locator(AutoMinorLocator(ny))
        ax.tick_params(which="both", top=False, right=False)

    # ---------------------------------------------------------- (a) result
    axA.axhline(LC_THEORY, color=BLUE, lw=1.6, zorder=2,
                label=r"$\rho(N)=1$")
    axA.axhline(LC_O1, color=CORAL, lw=1.6, ls=(0, (4, 2.4)), zorder=2,
                label=r"$C(3,2\lambda_c)=1$")
    axA.errorbar(opair, lcs, yerr=ses, fmt="o", ms=5.2, mfc=TEAL, mec="white",
                 mew=0.8, ecolor=TEAL, elinewidth=1.1, capsize=2.4,
                 capthick=0.9, zorder=4)
    axA.plot(opair, lcs, "-", color=TEAL, lw=1.3, zorder=3)
    axA.set_xlim(-0.06, 1.06)
    axA.set_xlabel(r"$o$")
    axA.set_ylabel(r"$\lambda_c$")
    minor(axA)
    axA.legend(loc="upper left", fontsize=8.4, handlelength=1.7,
               handletextpad=0.6, labelspacing=0.45, borderaxespad=0.7)
    tag(axA, "a")

    # ----------------------------------------------------- (b) measurement
    axB.axhline(0, color=MUTED, lw=0.8, zorder=1)
    for i, o in enumerate(os_):
        col = SEQ(NORM(o))
        lams, y, ye = lams_all[i], y_all[i], ye_all[i]
        a, b, _ = _wls(lams, y, ye)
        xi = -a / b
        xs = np.linspace(lams[0], xi, 40)
        axB.plot(xs, a + b * xs, "-", color=col, lw=1.4, zorder=3)
        axB.errorbar(lams, y, yerr=ye, fmt="o", ms=3.6, mfc=col, mec="white",
                     mew=0.6, ecolor=col, elinewidth=0.9, capsize=1.8,
                     capthick=0.7, zorder=4)
        axB.plot([xi], [0.0], "o", ms=4.6, mfc="white", mec=col, mew=1.1, zorder=5)
    axB.set_xlabel(r"$\lambda$")
    axB.set_ylabel(r"$1/\chi$")
    axB.set_ylim(bottom=-0.006)
    minor(axB)
    cax = axB.inset_axes([0.90, 0.56, 0.032, 0.36])
    cb = fig.colorbar(ScalarMappable(norm=NORM, cmap=SEQ), cax=cax)
    cb.set_ticks([0, 1])
    cb.ax.tick_params(labelsize=7.5, length=2.2, color=SEC)
    cb.outline.set_linewidth(0.5); cb.outline.set_edgecolor("#c6c5bf")
    cb.set_label(r"$o$", fontsize=8.6, labelpad=1.0)
    tag(axB, "b")

    # ------------------------------------------------------ (c) mechanism
    lam = np.linspace(0.0, 0.55, 260)
    th = np.array([3 * cascade_C(3, L, 1) for L in lam])        # what (2.7) counts
    ac = np.array([cascade_C(3, 2 * L, 1) for L in lam])        # the o=1 process
    axC.axhline(1.0, color=MUTED, lw=0.9, ls=(0, (4, 2.4)), zorder=2)
    axC.plot(lam, th, "-", color=BLUE, lw=1.8, zorder=3,
             label=r"$3\,C(3,\lambda)$")
    axC.plot(lam, ac, "-", color=CORAL, lw=1.8, zorder=3,
             label=r"$C(3,2\lambda)$")
    axC.plot([LC_THEORY], [1.0], "o", ms=5.4, mfc="white", mec=BLUE, mew=1.2, zorder=5)
    axC.plot([LC_O1], [1.0], "o", ms=5.4, mfc="white", mec=CORAL, mew=1.2, zorder=5)
    axC.set_xlim(0, 0.55); axC.set_ylim(0, 2.6)
    axC.set_xlabel(r"$\lambda$")
    axC.set_ylabel(r"$X\!\cdot\!C$")
    minor(axC)
    axC.legend(loc="upper left", fontsize=8.6, handlelength=1.7,
               handletextpad=0.6, labelspacing=0.5, borderaxespad=0.7)
    tag(axC, "c")

    fig.savefig("figure7_overlap.pdf", bbox_inches="tight")
    fig.savefig("figure7_overlap.png", bbox_inches="tight", dpi=210)
    print("saved figure7_overlap.png/.pdf")


if __name__ == "__main__":
    main()
