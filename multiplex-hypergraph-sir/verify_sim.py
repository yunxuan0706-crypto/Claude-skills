"""
Simulation validation for Figure 3.

Figure 3 compares two closures (mean field vs group closure), so on its own it
is theory against theory. This script supplies the missing ground truth in two
parts, chosen to match a structural fact: the sign-flip region always sits at
2.5-4.7 times lambda_c, i.e. it is *always supercritical*, so it cannot be
reached by a subcritical outbreak-size measurement. Therefore

  (1) the full next-generation bookkeeping N -- transmissibility T, the excess
      subtraction, and the cascade -- is validated where it can be, by exact
      Gillespie simulation on a real configuration-model multiplex hypergraph
      (with loops and finite N) in the subcritical regime, across m; and

  (2) the cascade factor C(m, lambda), which is the term that drives the flip,
      is validated by direct single-group Gillespie sampling *at the flip
      region's own parameters*.

Together these validate every ingredient the flip is built from.
"""
import numpy as np

from theory import lambda_c_rho, cascade_C, next_gen_matrix, degree_moments
from meanfield import ngm_switched, factors
from simulate import mean_outbreak_size
from recheck_mc import cascade_C_mc


def chi_from_matrix(N, g):
    """Mean outbreak size 1 + 1^T (I - N^T)^{-1} g for a seed's first
    generation g (branching prediction)."""
    M = N.shape[0]
    return 1.0 + np.linalg.solve(np.eye(M) - N.T, g).sum()


def chi_exact(P, m, lam):
    M = len(m)
    mean_k, _ = degree_moments(P, M)
    N = next_gen_matrix(P, m, lam)
    g = np.array([mean_k[b] * cascade_C(m[b], lam, 1) for b in range(M)])
    return chi_from_matrix(N, g)


def chi_meanfield(P, m, lam):
    """Same construction, but with the mean-field matrix and its
    mean-field first generation (m-1)*lam in place of the cascade."""
    M = len(m)
    mean_k, _ = degree_moments(P, M)
    N = ngm_switched(P, m, lam, s_T=0, s_D=0, s_C=0)
    g = np.array([mean_k[b] * (m[b] - 1) * lam for b in range(M)])
    return chi_from_matrix(N, g)


def ok(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name}  {detail}")
    return cond


if __name__ == "__main__":
    import json
    allpass = True
    P = {(3, 3): 0.5, (5, 5): 0.5}

    # ---------------------------------------------------------------- part 1
    print("=== 1. Gillespie on a real multiplex hypergraph vs the exact N ===")
    rows = []
    for m in (3, 5, 8, 12):
        lc = lambda_c_rho(P, (m, m))
        for frac in (0.6, 0.8):
            lam = frac * lc
            sim = mean_outbreak_size(P, (m, m), lam, N=40000, n_seeds=3000,
                                     n_graphs=3, seed=100 + m)
            ex = chi_exact(P, (m, m), lam)
            mf = chi_meanfield(P, (m, m), lam)
            z = (sim["chi"] - ex) / sim["sem"]
            rows.append(dict(m=m, frac=frac, lam=lam, chi_sim=sim["chi"],
                             sem=sim["sem"], chi_exact=ex, chi_mf=mf, z=z))
            print(f"  m={m:2d} lam={frac}lc : sim={sim['chi']:7.3f}+-{sim['sem']:.3f}"
                  f"  exact={ex:7.3f} ({z:+5.2f}s)  meanfield={mf:7.3f}"
                  f" ({(sim['chi']-mf)/sim['sem']:+6.1f}s)")
            allpass &= ok(f"   sim agrees with exact N (m={m}, {frac}lc)",
                          abs(z) < 3.5)

    n_mf_off = sum(1 for r in rows if abs((r["chi_sim"] - r["chi_mf"]) / r["sem"]) > 3.5)
    allpass &= ok("mean field is excluded by simulation at most points",
                  n_mf_off >= len(rows) - 2, f"{n_mf_off}/{len(rows)} points")

    # ---------------------------------------------------------------- part 2
    print("\n=== 2. Cascade C sampled directly AT the flip-region parameters ===")
    flip_pts = [(8, 0.15), (8, 0.30), (10, 0.20), (12, 0.30), (16, 0.15)]
    for m, lam in flip_pts:
        c_th = cascade_C(m, lam, 1)
        c_mc, se = cascade_C_mc(m, lam, 1, nruns=200_000)
        nz = abs(c_mc - c_th) / se
        fT, fD, fC, net = factors(P, (m, m), lam)
        allpass &= ok(f"C({m},{lam}) recursion == Gillespie", nz < 3.5,
                      f"C={c_th:.4f} MC={c_mc:.4f}+-{se:.4f} ({nz:.2f}s) "
                      f"| here net={net:.3f} (>1 => MF underestimates)")
        allpass &= ok(f"   this point is inside the flip region", net > 1.0)

    # ---------------------------------------------------------------- part 3
    print("\n=== 3. the flip region is always supercritical (structural) ===")
    worst = 1e9
    for Pd in [P, {(10, 10): 1.0}, {(2, 2): .9, (20, 20): .1}, {(3, 3): 1.0}]:
        for m in (8, 10, 12, 16):
            lc = lambda_c_rho(Pd, (m, m))
            ls = np.geomspace(1e-3, 20, 1500)
            nets = np.array([factors(Pd, (m, m), l)[3] for l in ls])
            w = ls[nets > 1]
            if len(w):
                worst = min(worst, w.min() / lc)
    allpass &= ok("flip window starts above lambda_c for every case",
                  worst > 1.0, f"closest approach = {worst:.2f} x lambda_c")

    json.dump(rows, open("figure3_sim.json", "w"), indent=2)
    print("\nwrote figure3_sim.json")
    print("\n=== ALL PASS ===" if allpass else "\n!!! SOME FAILURES !!!")
