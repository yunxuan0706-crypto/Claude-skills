"""Figure 2: time evolution of multiplex-hypergraph SIR -- group closure vs
exact Gillespie simulation (no mean field).

The group closure (theory.py) is integrated for S(t); one extra quadrature
Rdot = mu I splits I(t) from R(t). The ground truth is a direct continuous-time
(Gillespie) simulation of the microscopic process on a configuration-model
multiplex hypergraph, started from an epsilon fraction infected and averaged
over independent realisations. A Fenwick tree over group rates makes each
transmission draw O(log G), so a supercritical run at N=6000 is cheap.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from theory import Closure, lambda_c_rho
from simulate import build_multiplex

# ---------------------------------------------------------------- parameters
P = {(2, 2): 0.5, (3, 3): 0.5}
m = (3, 4)
N = 6000
EPS = 0.02
NRUNS = 400
LC = lambda_c_rho(P, m)
LAM = 1.6 * LC
TGRID = np.linspace(0.0, 20.0, 400)


# ---------------------------------------------------------------- simulation
def gillespie_traj(P, m, N, lam, eps, tgrid, rng, w=None, theta=None):
    """Exact continuous-time SIR trajectory; returns S(t), I(t) on tgrid."""
    M = len(m)
    w = np.ones(M) if w is None else np.asarray(w, float)
    theta = [1] * M if theta is None else theta
    beta = [lam * w[a] for a in range(M)]
    groups, _ = build_multiplex(P, m, N, rng)
    gmembers, glayer = [], []
    for a in range(M):
        for row in groups[a]:
            gmembers.append(np.asarray(row, np.int64))
            glayer.append(a)
    G = len(gmembers)
    glayer = np.array(glayer)
    membership = [[] for _ in range(N)]
    for gi, row in enumerate(gmembers):
        for v in row:
            membership[int(v)].append(gi)
    membership = [np.array(x, np.int64) for x in membership]

    st = np.zeros(N, np.int8)                      # 0=S 1=I 2=R
    st[rng.random(N) < eps] = 1
    n_inf = np.zeros(G, np.int64)
    n_sus = np.zeros(G, np.int64)
    for gi, row in enumerate(gmembers):
        s = st[row]
        n_inf[gi] = (s == 1).sum()
        n_sus[gi] = (s == 0).sum()

    def grate(gi):
        a = glayer[gi]
        return beta[a] * n_sus[gi] if n_inf[gi] >= theta[a] else 0.0

    bit = np.zeros(G + 1)

    def bupd(i, delta):
        i += 1
        while i <= G:
            bit[i] += delta
            i += i & (-i)

    def bprefix(i):
        s = 0.0
        while i > 0:
            s += bit[i]
            i -= i & (-i)
        return s

    def bfind(x):                                  # 0-based group with prefix>x
        pos, r = 0, x
        for k in range(G.bit_length(), -1, -1):
            nxt = pos + (1 << k)
            if nxt <= G and bit[nxt] <= r:
                pos = nxt
                r -= bit[nxt]
        return pos

    rate = np.array([grate(gi) for gi in range(G)])
    for gi in range(G):
        if rate[gi] > 0:
            bupd(gi, rate[gi])
    infected = set(np.where(st == 1)[0].tolist())

    def set_rate(gi):
        new = grate(gi)
        if new != rate[gi]:
            bupd(gi, new - rate[gi])
            rate[gi] = new

    S, I = np.empty(len(tgrid)), np.empty(len(tgrid))
    t, gt = 0.0, 0

    def record_until(tnow):
        nonlocal gt
        while gt < len(tgrid) and tgrid[gt] <= tnow:
            S[gt] = (st == 0).sum() / N
            I[gt] = len(infected) / N
            gt += 1

    record_until(0.0)
    while infected and gt < len(tgrid):
        rec = float(len(infected))
        tr = bprefix(G)
        tot = rec + tr
        if tot <= 0:
            break
        t += rng.exponential(1.0 / tot)
        record_until(t)
        if gt >= len(tgrid):
            break
        if rng.random() * tot < rec:                       # recovery
            arr = np.fromiter(infected, np.int64, len(infected))
            v = int(arr[rng.integers(len(arr))])
            infected.discard(v)
            st[v] = 2
            for gi in membership[v]:
                n_inf[gi] -= 1
                set_rate(gi)
        else:                                              # transmission
            gi = bfind(rng.random() * tr)
            row = gmembers[gi]
            sus = [int(x) for x in row if st[x] == 0]
            wnode = sus[rng.integers(len(sus))]
            st[wnode] = 1
            infected.add(wnode)
            for gj in membership[wnode]:
                n_inf[gj] += 1
                n_sus[gj] -= 1
                set_rate(gj)
    record_until(np.inf)
    return S, I


# ---------------------------------------------------------------- closure S,I
def closure_SI(tgrid):
    clo = Closure(P, m, LAM, eps=EPS)
    K = clo.dim

    def rhs(t, y):                                # x-compartments + R (Rdot=muI)
        x, R = y[:K], y[K]
        Phi, _ = clo.phi(x)
        S = (1 - EPS) * clo.gf.Psi(Phi)
        return np.concatenate([clo.rhs(t, x), [1 - S - R]])

    sol = solve_ivp(rhs, [0, tgrid[-1]],
                    np.concatenate([clo.initial(), [0.0]]),
                    method="LSODA", rtol=1e-9, atol=1e-12, dense_output=True)
    Y = sol.sol(tgrid)
    S = np.array([(1 - EPS) * clo.gf.Psi(clo.phi(Y[:K, j])[0])
                  for j in range(len(tgrid))])
    return S, 1 - S - Y[K, :]


def main():
    Sc, Ic = closure_SI(TGRID)
    rng = np.random.default_rng(20240601)
    Ss = np.empty((NRUNS, len(TGRID)))
    Is = np.empty((NRUNS, len(TGRID)))
    for r in range(NRUNS):
        Ss[r], Is[r] = gillespie_traj(P, m, N, LAM, EPS, TGRID, rng)
    Sm, Im = Ss.mean(0), Is.mean(0)
    Sci = 1.96 * Ss.std(0, ddof=1) / np.sqrt(NRUNS)
    Ici = 1.96 * Is.std(0, ddof=1) / np.sqrt(NRUNS)
    devS, devI = np.abs(Sm - Sc).max(), np.abs(Im - Ic).max()
    peak_sim, peak_clo = Im.max(), Ic.max()

    # ---- plot ----
    INK, SEC = "#20201e", "#565550"
    SIM, CLO = "#20201e", "#2f5fd0"
    BAND = "#c9c9c4"
    mpl.rcParams.update({
        "figure.dpi": 140, "savefig.dpi": 340,
        "font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 9.5,
        "axes.edgecolor": "#c6c5bf", "axes.linewidth": 0.7, "axes.labelcolor": INK,
        "xtick.color": SEC, "ytick.color": SEC, "text.color": INK,
        "axes.grid": True, "grid.color": "#ecebe5", "grid.linewidth": 0.55,
        "axes.axisbelow": True, "legend.frameon": False,
    })
    fig, (axS, axI) = plt.subplots(1, 2, figsize=(9.6, 3.35))
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.155, top=0.9, wspace=0.24)
    for ax, (ym, yci, yc, lab) in ((axS, (Sm, Sci, Sc, r"$S(t)$")),
                                   (axI, (Im, Ici, Ic, r"$I(t)$"))):
        ax.fill_between(TGRID, ym - yci, ym + yci, color=BAND, lw=0, zorder=1)
        ax.plot(TGRID, ym, "-", color=SIM, lw=1.5, zorder=3, label="simulation")
        ax.plot(TGRID, yc, "--", color=CLO, lw=1.6, dashes=(4, 2.4), zorder=4,
                label="group closure")
        ax.set_xlim(0, 20)
        ax.set_xlabel(r"time $t$")
        ax.set_ylabel(lab)
        ax.tick_params(which="both", top=False, right=False)
    axS.legend(loc="center right", fontsize=8.6, handlelength=1.8,
               handletextpad=0.6, labelspacing=0.5)
    axS.text(-0.16, 1.04, "a", transform=axS.transAxes, fontsize=13,
             fontweight="bold", va="bottom")
    axI.text(-0.16, 1.04, "b", transform=axI.transAxes, fontsize=13,
             fontweight="bold", va="bottom")
    fig.savefig("figure2_evolution.pdf", bbox_inches="tight")
    fig.savefig("figure2_evolution.png", bbox_inches="tight", dpi=210)
    print(f"saved figure2_evolution.png/.pdf  (N={N}, {NRUNS} runs, lam=1.6 lc)")
    print(f"max |S_sim-S_clo| = {devS:.4f}   max |I_sim-I_clo| = {devI:.4f}")
    print(f"peak I: sim {peak_sim:.4f}  closure {peak_clo:.4f}  "
          f"at t={TGRID[Im.argmax()]:.2f}/{TGRID[Ic.argmax()]:.2f}")


if __name__ == "__main__":
    main()
