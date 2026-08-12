"""
Independent data checks for the subcritical extrapolation.

Key check: the closure-ODE mean outbreak size  chi_ODE(lambda) = lim_{eps->0} R(inf)/eps
must equal the multi-type branching prediction
    chi_branch = 1 + 1^T (I - N^T)^{-1} g ,   g_b = <k_b> C(m_b, lam_b, theta_b),
which is the outline's  1^T (I-N)^{-1} v0  written for a uniformly random seed.
These share NO code with each other (ODE integration vs. a matrix inverse), so
agreement validates the final-state data the figure is built on.
"""
import numpy as np
from theory import (Closure, next_gen_matrix, degree_moments, cascade_C,
                    lambda_c_rho, GenFun)
from lambda_c_extrap import chi_subcritical, CONFIGS


def chi_branch(P, m, lam, w=None, theta=None):
    M = len(m)
    if w is None: w = np.ones(M)
    if theta is None: theta = np.ones(M, dtype=int)
    mean_k, _ = degree_moments(P, M)
    N = next_gen_matrix(P, m, lam, w, theta)
    Cb = np.array([cascade_C(m[b], lam*w[b], theta[b]) for b in range(M)])
    g = mean_k * Cb                                   # seed's first generation
    x = np.linalg.solve(np.eye(M) - N.T, g)           # (I-N^T)^{-1} g
    return 1.0 + x.sum()


print("=== Check 1: closure-ODE outbreak size  vs  branching formula ===")
worst = 0.0
for name, P, m, theta in CONFIGS:
    lc = lambda_c_rho(P, m, theta=theta)
    for f in (0.80, 0.90, 0.95):
        lam = f*lc
        chi_o, _ = chi_subcritical(P, m, lam, theta=theta, eps=1e-4)
        chi_b = chi_branch(P, m, lam, theta=theta)
        rel = abs(chi_o - chi_b)/chi_b
        worst = max(worst, rel)
    print(f"  {name:16s} lam=0.95lc: chi_ODE={chi_o:8.3f}  "
          f"chi_branch={chi_b:8.3f}  rel.diff={rel:.2e}")
print(f"  worst relative difference over all configs/lambdas: {worst:.2e}\n")

print("=== Check 2: eps->0 Richardson stability (halving the seed again) ===")
for name, P, m, theta in CONFIGS[:3]:
    lc = lambda_c_rho(P, m, theta=theta)
    lam = 0.95*lc
    c1, _ = chi_subcritical(P, m, lam, theta=theta, eps=1e-4)
    c2, _ = chi_subcritical(P, m, lam, theta=theta, eps=5e-5)
    print(f"  {name:16s} chi(eps=1e-4)={c1:.5f}  chi(eps=5e-5)={c2:.5f}  "
          f"rel={abs(c1-c2)/c2:.1e}")

print("\n=== Check 3: final-state convergence in t_max ===")
for name, P, m, theta in CONFIGS[:3]:
    lc = lambda_c_rho(P, m, theta=theta)
    lam = 0.995*lc            # closest-to-threshold point (slowest)
    for tmax in (2e4, 6e4, 1.2e5):
        c, leftover = chi_subcritical(P, m, lam, theta=theta, eps=1e-4, tmax=tmax)
        print(f"  {name:16s} lam=0.995lc tmax={tmax:.0e}: chi={c:.4f} "
              f"leftover={leftover:.1e}")

print("\n=== Check 4: figure_data.json sanity (monotone, positive, fit) ===")
import json
d = json.load(open("figure_data.json"))
allok = True
for r in d:
    ys = np.array(r["ys"]); lams = np.array(r["lams"])
    mono = np.all(np.diff(ys) < 0)               # 1/chi decreasing toward lc
    pos = np.all(ys > 0)
    a, b = r["fit"]
    x0 = -a/b
    close = abs(x0 - r["lc_extrap"]) < 1e-9
    ok = mono and pos and close
    allok &= ok
    print(f"  {r['name']:16s} decreasing={mono} positive={pos} "
          f"fit_consistent={close}")
print("  ALL SANE" if allok else "  PROBLEM")
