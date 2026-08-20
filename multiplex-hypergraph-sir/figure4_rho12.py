"""Figure 4: how inter-layer participation correlation moves the threshold.

A node's participation is its layer-degree vector (k^(1), k^(2)) -- how many
groups it joins in each layer. We hold the marginals fixed (each layer degree
is 2 or 4 with equal probability, mean 3) and vary only the Pearson correlation
rho12 between the two layer degrees:

    P(2,2) = P(4,4) = (1+rho12)/4 ,   P(2,4) = P(4,2) = (1-rho12)/4 .

(With degrees {2,4} the per-layer variance is 1, so the free parameter equals
rho12 exactly.) Positive rho12 means the same nodes are hubs in both layers,
which concentrates the cross-layer spreading capacity (the off-diagonal
<k^1 k^2> in the next-generation matrix) and LOWERS the threshold; anti-
correlation raises it.

Three panels tell the full story:
  (a) result      -- exact threshold rho(N)=1 vs a bootstrap of the Gillespie
                     measurement of lambda_c at six correlations;
  (b) measurement -- the subcritical 1/chi -> 0 extrapolation behind each box,
                     coloured by rho12 (its x-intercept IS lambda_c);
  (c) mechanism   -- the next-generation spectral radius rho(N) at fixed lambda
                     rising through unity as rho12 grows, i.e. why lambda_c falls.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import AutoMinorLocator
from theory import lambda_c_rho, rho_of_lambda

M_GROUPS = (3, 3)
N_SIM = 20000
NSEEDS = 5000
NGRAPHS = 4
WINDOW = np.array([0.72, 0.77, 0.82, 0.87, 0.92])   # fractions of lambda_c
RHO_SIM = [-1.0, -0.6, -0.2, 0.2, 0.6, 1.0]

# ---------------------------------------------------------------- palette
INK, SEC, MUTED = "#20201e", "#565550", "#9b9a93"
TEAL, CORAL, BLUE = "#1C9B8E", "#E76F51", "#2F5FD0"      # shared jewel triad
DIV = LinearSegmentedColormap.from_list("rho12", [TEAL, "#eceae3", CORAL])
NORM = Normalize(-1.0, 1.0)


def Pfam(rho):
    """Joint layer-degree distribution with fixed {2,4} marginals, corr = rho."""
    tab = [((2, 2), (1 + rho) / 4), ((4, 4), (1 + rho) / 4),
           ((2, 4), (1 - rho) / 4), ((4, 2), (1 - rho) / 4)]
    return {k: p for k, p in tab if p > 1e-12}


def _wls(lams, y, ye):
    """Weighted least-squares intercept a and slope b of y ~ a + b*lam."""
    W = 1.0 / ye ** 2
    X = np.vstack([np.ones_like(lams), lams]).T
    cov = np.linalg.inv(X.T @ (W[:, None] * X))
    a, b = cov @ (X.T @ (W * y))
    return a, b, cov


def _xintercept(lams, y, ye):
    a, b, cov = _wls(lams, y, ye)
    g = np.array([-1.0 / b, a / b ** 2])             # d(-a/b)
    return -a / b, float(np.sqrt(g @ cov @ g))


def lambda_c_sim(rho, seed=0):
    """Simulation lambda_c by 1/chi -> 0 extrapolation.
    Returns (lambda_c, se, lams, y=1/chi, ye)."""
    from simulate import mean_outbreak_size
    P = Pfam(rho)
    lams = WINDOW * lambda_c_rho(P, M_GROUPS)
    y, ye = [], []
    for L in lams:
        d = mean_outbreak_size(P, M_GROUPS, L, N=N_SIM, n_seeds=NSEEDS,
                               n_graphs=NGRAPHS, seed=seed)
        y.append(1.0 / d["chi"])
        ye.append(d["sem"] / d["chi"] ** 2)          # delta method on 1/chi
    y, ye = np.array(y), np.array(ye)
    lcs, se = _xintercept(lams, y, ye)
    return lcs, se, lams, y, ye


def bootstrap_lc(lams, y, ye, B=4000, rng=None):
    """Parametric bootstrap of the fitted threshold: resample each 1/chi
    ~ Normal(mean, se) and refit, returning the B lambda_c draws."""
    rng = np.random.default_rng(0) if rng is None else rng
    out = np.empty(B)
    for j in range(B):
        out[j] = _xintercept(lams, y + ye * rng.standard_normal(len(y)), ye)[0]
    return out


def main():
    import json, os
    rhos = np.linspace(-1, 1, 41)
    lc_theory = np.array([lambda_c_rho(Pfam(r), M_GROUPS) for r in rhos])

    # ---- simulation, per-lambda data cached in figure4_data.json ----
    if os.path.exists("figure4_data.json"):
        d = json.load(open("figure4_data.json"))
        xs, ys, es = (np.array(d["rho_sim"]), np.array(d["lc_sim"]),
                      np.array(d["se_sim"]))
        lams_all = [np.array(v) for v in d["lams_all"]]
        y_all = [np.array(v) for v in d["y_all"]]
        ye_all = [np.array(v) for v in d["ye_all"]]
    else:
        xs, ys, es, lams_all, y_all, ye_all = [], [], [], [], [], []
        for i, r in enumerate(RHO_SIM):
            lcs, se, lams, y, ye = lambda_c_sim(r, seed=100 + i)
            xs.append(r); ys.append(lcs); es.append(se)
            lams_all.append(lams); y_all.append(y); ye_all.append(ye)
            print(f"rho12={r:+.2f}: theory={lambda_c_rho(Pfam(r), M_GROUPS):.5f}"
                  f"  sim={lcs:.5f}+-{se:.5f}")
        xs, ys, es = np.array(xs), np.array(ys), np.array(es)
        json.dump({"rho_sim": xs.tolist(), "lc_sim": ys.tolist(),
                   "se_sim": es.tolist(),
                   "lams_all": [v.tolist() for v in lams_all],
                   "y_all": [v.tolist() for v in y_all],
                   "ye_all": [v.tolist() for v in ye_all],
                   "N_sim": N_SIM, "n_seeds": NSEEDS, "window": WINDOW.tolist()},
                  open("figure4_data.json", "w"), indent=2)

    # ---- bootstrap the lambda_c distribution at each rho12 ----
    rng = np.random.default_rng(2024)
    boots = [bootstrap_lc(lams_all[i], y_all[i], ye_all[i], B=4000, rng=rng)
             for i in range(len(xs))]

    # ================================================================= style
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
    fig.subplots_adjust(left=0.058, right=0.985, bottom=0.15, top=0.905, wspace=0.30)

    def tag(ax, s):
        ax.text(-0.20, 1.045, s, transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="bottom", ha="left", color=INK)

    def minorticks(ax, nx=2, ny=2):
        ax.xaxis.set_minor_locator(AutoMinorLocator(nx))
        ax.yaxis.set_minor_locator(AutoMinorLocator(ny))
        ax.tick_params(which="both", top=False, right=False)

    # ============================================================ (a) result
    MED = "#a23c1a"      # dark coral for box edges/median, reads as coral family
    axA.plot(rhos, lc_theory, "-", color=BLUE, lw=1.9, zorder=2,
             label=r"theory, $\rho(N)=1$")
    bp = axA.boxplot(boots, positions=xs, widths=0.115, patch_artist=True,
                     whis=(5, 95), showfliers=False, showcaps=True,
                     manage_ticks=False, zorder=4,
                     medianprops=dict(color=MED, lw=1.4),
                     boxprops=dict(lw=0.9, edgecolor=MED),
                     whiskerprops=dict(color=MED, lw=0.9),
                     capprops=dict(color=MED, lw=0.9))
    for patch in bp["boxes"]:
        patch.set_facecolor("#f4c9b6")     # soft coral fill
        patch.set_alpha(0.9)
    axA.set_xlim(-1.15, 1.15)
    axA.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
    axA.set_xlabel(r"inter-layer participation correlation $\rho_{12}$")
    axA.set_ylabel(r"epidemic threshold $\lambda_c$")
    minorticks(axA)
    from matplotlib.patches import Patch
    box_proxy = Patch(facecolor="#f4c9b6", edgecolor=MED, lw=0.9,
                      label="simulation (bootstrap)")
    axA.legend(handles=[axA.lines[0], box_proxy], loc="upper right", fontsize=8.6,
               handlelength=1.7, handletextpad=0.6, labelspacing=0.55)
    axA.text(0.035, 0.06, "correlated participation\nlowers the threshold",
             transform=axA.transAxes, fontsize=8.4, color=SEC, ha="left",
             va="bottom", linespacing=1.4)
    tag(axA, "a")

    # ==================================================== (b) the measurement
    axB.axhline(0, color=MUTED, lw=0.8, zorder=1)
    for i, r in enumerate(RHO_SIM):
        col = DIV(NORM(r))
        lams, y, ye = lams_all[i], y_all[i], ye_all[i]
        a, b, _ = _wls(lams, y, ye)
        xi = -a / b
        xln = np.linspace(lams[0], xi, 40)
        axB.plot(xln, a + b * xln, "-", color=col, lw=1.4, zorder=3)
        axB.errorbar(lams, y, yerr=ye, fmt="o", ms=3.8, mfc=col, mec="white",
                     mew=0.6, ecolor=col, elinewidth=0.9, capsize=1.8,
                     capthick=0.7, zorder=4)
        axB.plot([xi], [0.0], "o", ms=4.6, mfc="white", mec=col, mew=1.1, zorder=5)
    axB.set_xlabel(r"transmissibility $\lambda$")
    axB.set_ylabel(r"$1/\chi$   (inverse mean cluster size)")
    axB.set_ylim(bottom=-0.006)
    minorticks(axB)
    # slim vertical colourbar (rho12) tucked into the empty lower-left wedge
    cax = axB.inset_axes([0.085, 0.12, 0.035, 0.42])
    cb = fig.colorbar(ScalarMappable(norm=NORM, cmap=DIV), cax=cax)
    cb.set_ticks([-1, 0, 1])
    cb.ax.tick_params(labelsize=7.5, length=2.2, color=SEC)
    cb.outline.set_linewidth(0.5); cb.outline.set_edgecolor("#c6c5bf")
    cb.set_label(r"$\rho_{12}$", fontsize=8.6, labelpad=1.0)
    axB.text(0.62, 0.9, r"$1/\chi\!\to\!0$ at $\lambda_c$", transform=axB.transAxes,
             fontsize=8.4, color=SEC, ha="left", va="top")
    tag(axB, "b")

    # ====================================================== (c) the mechanism
    lc0 = lambda_c_rho(Pfam(0.0), M_GROUPS)
    shades = ["#9db8e0", BLUE, "#1c3a80"]
    fracs = [0.92, 1.0, 1.08]
    for frac, sh in zip(fracs, shades):
        lam = frac * lc0
        rr = np.array([rho_of_lambda(lam, Pfam(r), M_GROUPS) for r in rhos])
        axC.plot(rhos, rr, "-", color=sh, lw=1.7, zorder=3,
                 label=rf"$\lambda={frac:.2f}\,\lambda_c(0)$")
    axC.axhline(1.0, color=CORAL, lw=1.2, ls=(0, (4, 2.4)), zorder=2)
    axC.plot([0.0], [1.0], "o", ms=5.2, mfc="white", mec=CORAL, mew=1.2, zorder=5)
    axC.text(0.98, 1.0, r"$\rho(N)=1$", color=CORAL, fontsize=8.6, ha="right",
             va="bottom", transform=axC.get_yaxis_transform())
    axC.set_xlim(-1.05, 1.05)
    axC.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
    axC.set_xlabel(r"inter-layer participation correlation $\rho_{12}$")
    axC.set_ylabel(r"spectral radius $\rho(N)$   (fixed $\lambda$)")
    minorticks(axC)
    axC.legend(loc="upper left", fontsize=8.2, handlelength=1.6,
               handletextpad=0.55, labelspacing=0.45, borderaxespad=0.6)
    tag(axC, "c")

    fig.savefig("figure4_rho12.pdf", bbox_inches="tight")
    fig.savefig("figure4_rho12.png", bbox_inches="tight", dpi=210)
    print("saved figure4_rho12.png/.pdf  (3-panel)")
    print(f"lambda_c: {lc_theory[0]:.4f} (rho12=-1) -> {lc_theory[-1]:.4f} "
          f"(rho12=+1), a {(1-lc_theory[-1]/lc_theory[0])*100:.1f}% drop")
    # sanity: the fixed-lambda=lc0 curve must cross unity at rho12=0
    print(f"rho(N; lambda=lc0) at rho12=0 : "
          f"{rho_of_lambda(lc0, Pfam(0.0), M_GROUPS):.6f}  (should be 1)")


if __name__ == "__main__":
    main()
