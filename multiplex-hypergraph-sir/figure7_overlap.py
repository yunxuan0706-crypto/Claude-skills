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
WINDOW = np.array([0.72, 0.77, 0.82, 0.87, 0.92])
O_SIM = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

P_FAM = {(KDEG, KDEG): 1.0}
LC_THEORY = lambda_c_rho(P_FAM, M_GROUPS)                       # 0.186141
LC_O1 = brentq(lambda L: cascade_C(3, L, 1) - 1.0, 1e-9, 60, xtol=1e-14) / 2


def chi_at(o, lam, n_seeds, seed, cap=None):
    """Mean subcritical outbreak size, averaged over NGRAPHS realisations."""
    rng = np.random.default_rng(seed)
    cap = cap or max(2000, N_SIM // 10)
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
    """Two passes: locate lambda_c coarsely, then measure on a window placed
    just below it."""
    guess = LC_THEORY + o * (LC_O1 - LC_THEORY)          # bracketing interpolation
    for npass, nseed in ((0, 1200), (1, 5000)):
        lams = WINDOW * guess
        y, ye = [], []
        for L in lams:
            c, se = chi_at(o, L, nseed, seed + 17 * npass + int(1000 * L))
            y.append(1.0 / c)
            ye.append(se / c ** 2)
        y, ye = np.array(y), np.array(ye)
        lc, selc = _xintercept(lams, y, ye)
        guess = lc                                        # re-centre for pass 2
    return lc, selc, lams, y, ye


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
            lc, se, lams, y, ye = lambda_c_of_o(o, seed=500 + 31 * i)
            op = realised_overlap(o, seed=9000 + i)
            os_.append(o); lcs.append(lc); ses.append(se); opair.append(op)
            lams_all.append(lams); y_all.append(y); ye_all.append(ye)
            print(f"o={o:.2f} (pair {op:.4f}): lambda_c={lc:.5f}+-{se:.5f}"
                  f"  theory={LC_THEORY:.5f}  ratio={lc/LC_THEORY:.3f}", flush=True)
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
    return os_, opair, lcs, ses, lams_all, y_all, ye_all


if __name__ == "__main__":
    main()
