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
correlation raises it. The exact threshold rho(N)=1 is compared against a
direct Gillespie measurement of lambda_c (subcritical mean-cluster-size
extrapolation, 1/chi -> 0), with fixed marginals throughout.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from theory import lambda_c_rho
from simulate import mean_outbreak_size

M_GROUPS = (3, 3)
N_SIM = 20000
NSEEDS = 5000
NGRAPHS = 4
WINDOW = np.array([0.72, 0.77, 0.82, 0.87, 0.92])   # fractions of lambda_c
RHO_SIM = [-1.0, -0.6, -0.2, 0.2, 0.6, 1.0]


def Pfam(rho):
    """Joint layer-degree distribution with fixed {2,4} marginals, corr = rho."""
    tab = [((2, 2), (1 + rho) / 4), ((4, 4), (1 + rho) / 4),
           ((2, 4), (1 - rho) / 4), ((4, 2), (1 - rho) / 4)]
    return {k: p for k, p in tab if p > 1e-12}


def _xintercept(lams, y, ye):
    W = 1.0 / ye ** 2
    X = np.vstack([np.ones_like(lams), lams]).T
    cov = np.linalg.inv(X.T @ (W[:, None] * X))
    a, b = cov @ (X.T @ (W * y))
    g = np.array([-1.0 / b, a / b ** 2])             # d(-a/b)
    return -a / b, float(np.sqrt(g @ cov @ g))


def lambda_c_sim(rho, seed=0):
    """Simulation lambda_c by 1/chi -> 0 extrapolation.
    Returns (lambda_c, se, lams, y=1/chi, ye)."""
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

    # ---- plot: exact curve + bootstrap box plots ----
    INK, SEC = "#20201e", "#565550"
    THEORY, BOX, MED = "#2f5fd0", "#8fb4f2", "#1a3a7a"
    mpl.rcParams.update({
        "figure.dpi": 140, "savefig.dpi": 340,
        "font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 10,
        "axes.edgecolor": "#c6c5bf", "axes.linewidth": 0.7, "axes.labelcolor": INK,
        "xtick.color": SEC, "ytick.color": SEC, "text.color": INK,
        "axes.grid": True, "grid.color": "#ecebe5", "grid.linewidth": 0.55,
        "axes.axisbelow": True, "legend.frameon": False,
    })
    fig, ax = plt.subplots(figsize=(5.6, 4.1))
    fig.subplots_adjust(left=0.15, right=0.965, bottom=0.135, top=0.955)
    ax.plot(rhos, lc_theory, "-", color=THEORY, lw=1.9, zorder=2,
            label=r"theory, $\rho(N)=1$")
    bp = ax.boxplot(boots, positions=xs, widths=0.115, patch_artist=True,
                    whis=(5, 95), showfliers=False, showcaps=True,
                    manage_ticks=False, zorder=4,
                    medianprops=dict(color=MED, lw=1.4),
                    boxprops=dict(lw=0.9, edgecolor=MED),
                    whiskerprops=dict(color=MED, lw=0.9),
                    capprops=dict(color=MED, lw=0.9))
    for patch in bp["boxes"]:
        patch.set_facecolor(BOX)
        patch.set_alpha(0.85)
    ax.set_xlim(-1.15, 1.15)
    ax.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
    ax.set_xlabel(r"inter-layer participation correlation $\rho_{12}$")
    ax.set_ylabel(r"epidemic threshold $\lambda_c$")
    ax.tick_params(which="both", top=False, right=False)
    from matplotlib.patches import Patch
    box_proxy = Patch(facecolor=BOX, edgecolor=MED, lw=0.9,
                      label="simulation (bootstrap)")
    ax.legend(handles=[ax.lines[0], box_proxy], loc="upper right", fontsize=9,
              handlelength=1.8, handletextpad=0.6, labelspacing=0.6)
    ax.text(0.04, 0.075, "correlated participation\nlowers the threshold",
            transform=ax.transAxes, fontsize=8.8, color=SEC, ha="left",
            va="bottom", linespacing=1.4)
    fig.savefig("figure4_rho12.pdf", bbox_inches="tight")
    fig.savefig("figure4_rho12.png", bbox_inches="tight", dpi=210)
    print("saved figure4_rho12.png/.pdf")
    print(f"lambda_c: {lc_theory[0]:.4f} (rho12=-1) -> {lc_theory[-1]:.4f} "
          f"(rho12=+1), a {(1-lc_theory[-1]/lc_theory[0])*100:.1f}% drop")


if __name__ == "__main__":
    main()
