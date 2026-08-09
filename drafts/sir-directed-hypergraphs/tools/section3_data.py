#!/usr/bin/env python3
"""Publication-grade re-run of the Section III measurements.

The first pass (ebcm_directed.py) established the physics but is not fit to
plot, for three reasons that this script fixes:

  1. No uncertainty on any simulated quantity. Every number was a bare mean
     over stochastic runs. Here every simulated point carries an error bar,
     taken across independent STRUCTURES where a structure ensemble is what
     the claim is about, and across runs otherwise.

  2. The finite-size claim was confounded. N=1000 and N=4000 used different
     bi-degree draws (kappa 1.910 vs 1.988), hence different lambda, and
     different run counts (600 vs 300). Any of the three could produce the
     observed drop. Here the generating parameters, the run count and the
     ratio lambda/lambda_c are all held fixed, and only N moves.

  3. The r_io scan used one realisation per point, so kappa carried unquantified
     sampling noise. Here each point averages several independent sequences.

Output is a JSON blob consumed by section3_figures.py, so the figures can be
restyled without re-running the simulations.
"""
import json
import math
import os
import random
import statistics as st
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ebcm_directed import (PGF, bidegree, config_model, degrees, ebcm,
                           gillespie_traj, meanfield)
from feasibility_check import index

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "data", "section3.json")


def sem(v):
    return math.sqrt(st.variance(v) / len(v)) if len(v) > 1 else 0.0


def one_structure(N, tau, eta, md, corr, seed):
    rng = random.Random(seed)
    kin, kout = bidegree(N, tau, eta, md, corr, rng)
    edges = config_model(kin, kout, tau, eta, rng)
    ki, ko = degrees(N, edges)
    return edges, ki, ko, PGF(ki, ko)


def measured_rio(ki, ko):
    mi, mo = st.mean(ki), st.mean(ko)
    si, so = st.pstdev(ki), st.pstdev(ko)
    if si == 0 or so == 0:
        return float("nan")
    return (sum(a * b for a, b in zip(ki, ko)) / len(ki) - mi * mo) / (si * so)


# ---------------------------------------------------------------- panel 1 ----
def finite_size(Ns=(500, 1000, 2000, 4000), reps=4, runs=300,
                tau=2, eta=2, md=2, ratio=1.6, eps=0.01, seed=1000):
    """EBCM and mean-field residuals against exact simulation, versus N.

    Everything except N is held fixed: same generating parameters, same number
    of Gillespie runs, same lambda/lambda_c. reps independent structures per N
    give the error bar.
    """
    grid = [0.5 * k for k in range(1, 25)]           # t = 0.5 .. 12
    out = []
    for N in Ns:
        eb, mf = [], []
        for r in range(reps):
            edges, ki, ko, p = one_structure(N, tau, eta, md, 0.0,
                                             seed + 131 * r + N)
            lc = 1.0 / (tau * p.kappa - 1)
            lam = ratio * lc
            tails, heads, tails_of = index(N, edges)
            acc = [0.0] * len(grid)
            for s in range(runs):
                rr = random.Random(7_000_000 + 9173 * s + 17 * r + N)
                seeds = [v for v in range(N) if rr.random() < eps]
                tr = gillespie_traj(N, tails, heads, tails_of, lam, 1, rr,
                                    seeds, grid)
                for j, (S, _I, _R) in enumerate(tr):
                    acc[j] += S
            sim = [a / runs for a in acc]
            e = ebcm(p, tau, lam, eps=eps, tmax=13.0, dt=0.005)
            m = meanfield(N, edges, ki, ko, tau, lam, eps=eps, tmax=13.0,
                          dt=0.005)
            pick = lambda sol, t: sol[min(range(len(sol)),
                                          key=lambda k: abs(sol[k][0] - t))][1]
            eb.append(st.mean(abs(pick(e, t) - sim[j]) for j, t in enumerate(grid)))
            mf.append(st.mean(abs(pick(m, t) - sim[j]) for j, t in enumerate(grid)))
        out.append({"N": N, "ebcm": st.mean(eb), "ebcm_sem": sem(eb),
                    "mf": st.mean(mf), "mf_sem": sem(mf), "reps": reps,
                    "runs": runs})
        print(f"  N={N:5d}: EBCM {out[-1]['ebcm']:.4f}+/-{out[-1]['ebcm_sem']:.4f}"
              f"   MF {out[-1]['mf']:.4f}+/-{out[-1]['mf_sem']:.4f}", flush=True)
    return out


# ---------------------------------------------------------------- panel 2 ----
def threshold_sweep(N=6000, tau=2, eta=2, md=2, runs=400, eps=0.005,
                    ratios=(0.6, 0.8, 0.9, 1.0, 1.1, 1.3, 1.6, 2.0, 2.4),
                    seed=2000):
    """Final size across the threshold, with a run-level error bar."""
    edges, ki, ko, p = one_structure(N, tau, eta, md, 0.0, seed)
    lc = 1.0 / (tau * p.kappa - 1)
    tails, heads, tails_of = index(N, edges)
    rows = []
    for r in ratios:
        lam = r * lc
        vals = []
        for s in range(runs):
            rr = random.Random(3_000_000 + 6151 * s + int(1000 * r))
            seeds = [v for v in range(N) if rr.random() < eps]
            vals.append(gillespie_traj(N, tails, heads, tails_of, lam, 1, rr,
                                       seeds, [400.0])[-1][2])
        e = ebcm(p, tau, lam, eps=eps, tmax=250.0, dt=0.01)[-1][3]
        m = meanfield(N, edges, ki, ko, tau, lam, eps=eps, tmax=250.0,
                      dt=0.01)[-1][3]
        rows.append({"ratio": r, "sim": st.mean(vals), "sim_sem": sem(vals),
                     "ebcm": e, "mf": m})
        print(f"  lam/lc={r:.2f}: sim {rows[-1]['sim']:.4f}"
              f"+/-{rows[-1]['sim_sem']:.4f}  ebcm {e:.4f}  mf {m:.4f}",
              flush=True)
    return {"N": N, "kappa": p.kappa, "lc": lc, "lc_mf": 1 / (tau * p.kappa),
            "runs": runs, "eps": eps, "rows": rows}


# ---------------------------------------------------------------- panel 3 ----
def rio_scan(N=3000, tau=2, eta=2, md=2, reps=5,
             corrs=(-0.9, -0.7, -0.5, -0.3, 0.0, 0.3, 0.5, 0.7, 0.9),
             seed=3000):
    """lambda_c versus the realised r_io, several sequences per target."""
    rows = []
    for c in corrs:
        rs, lcs, kap = [], [], []
        mom_ki, mom_si, mom_so, mom_nm, lc38 = [], [], [], [], []
        for r in range(reps):
            e_, ki, ko, p = one_structure(N, tau, eta, md, c,
                                          seed + 977 * r + int(100 * c))
            if tau * p.kappa <= 1:
                continue
            rio = measured_rio(ki, ko)
            rs.append(rio)
            kap.append(p.kappa)
            lcs.append(1.0 / (tau * p.kappa - 1))
            si, so = st.pstdev(ki), st.pstdev(ko)
            nm = N / len(e_)
            mom_ki.append(st.mean(ki)); mom_si.append(si)
            mom_so.append(so); mom_nm.append(nm)
            lc38.append(1.0 / (tau * st.mean(ki) + nm * rio * si * so - 1))
        if not lcs:
            continue
        rows.append({"target": c, "rio": st.mean(rs), "rio_sem": sem(rs),
                     "kappa": st.mean(kap), "lc": st.mean(lcs),
                     "lc_sem": sem(lcs), "reps": len(lcs),
                     "ki_mean": st.mean(mom_ki), "si": st.mean(mom_si),
                     "so": st.mean(mom_so), "NoverM": st.mean(mom_nm),
                     # lambda_c re-expressed through eq (38); equal to the value
                     # above only if kappa = <k_in> + r_io*s_i*s_o/<k_out>, so
                     # this is an algebraic check of that rewriting, not an
                     # independent measurement of the threshold
                     "lc_eq38": st.mean(lc38)})
        print(f"  r_io={rows[-1]['rio']:+.3f}+/-{rows[-1]['rio_sem']:.3f}"
              f"  kappa={rows[-1]['kappa']:.3f}  lc={rows[-1]['lc']:.4f}"
              f"+/-{rows[-1]['lc_sem']:.4f}", flush=True)
    return rows


# ---------------------------------------------------------------- panel 4 ----
def trajectory(N=6000, tau=2, eta=2, md=2, runs=400, ratio=1.6, eps=0.01,
               seed=4000):
    """S(t) and I(t): simulation with a run-level band, versus both closures."""
    edges, ki, ko, p = one_structure(N, tau, eta, md, 0.0, seed)
    lc = 1.0 / (tau * p.kappa - 1)
    lam = ratio * lc
    grid = [0.25 * k for k in range(0, 57)]
    tails, heads, tails_of = index(N, edges)
    accS = [[] for _ in grid]
    accI = [[] for _ in grid]
    for s in range(runs):
        rr = random.Random(5_000_000 + 4441 * s)
        seeds = [v for v in range(N) if rr.random() < eps]
        tr = gillespie_traj(N, tails, heads, tails_of, lam, 1, rr, seeds, grid)
        for j, (S, I, _R) in enumerate(tr):
            accS[j].append(S)
            accI[j].append(I)
    e = ebcm(p, tau, lam, eps=eps, tmax=14.5, dt=0.005)
    m = meanfield(N, edges, ki, ko, tau, lam, eps=eps, tmax=14.5, dt=0.005)
    pick = lambda sol, t, i: sol[min(range(len(sol)),
                                     key=lambda k: abs(sol[k][0] - t))][i]
    return {"N": N, "lam": lam, "lc": lc, "ratio": ratio, "runs": runs,
            "t": grid,
            "S_sim": [st.mean(v) for v in accS],
            "S_sem": [sem(v) for v in accS],
            "I_sim": [st.mean(v) for v in accI],
            "I_sem": [sem(v) for v in accI],
            "S_ebcm": [pick(e, t, 1) for t in grid],
            "I_ebcm": [pick(e, t, 2) for t in grid],
            "S_mf": [pick(m, t, 1) for t in grid],
            "I_mf": [pick(m, t, 2) for t in grid]}


# ------------------------------------------------------ deterministic tables --
def convergence_tables():
    """Internal-consistency data: deterministic, so no error bars apply."""
    from ebcm_directed import states
    rng = random.Random(23)
    _e, ki, ko, p = one_structure(1500, 2, 2, 2, 0.0, 23)
    pred = 1.0 / (2 * p.kappa - 1)
    bis = []
    for tmax in (200.0, 400.0, 1600.0, 6400.0):
        lo, hi = 1e-3, 50.0
        while hi - lo > 1e-5 * hi:
            mid = math.sqrt(lo * hi)
            if ebcm(p, 2, mid, eps=1e-12, tmax=tmax, dt=0.05)[-1][3] > 1e-8:
                hi = mid
            else:
                lo = mid
        meas = math.sqrt(lo * hi)
        bis.append({"tmax": tmax, "measured": meas, "predicted": pred,
                    "rel": abs(meas - pred) / pred})
        print(f"  tmax={tmax:6.0f}: rel.err {bis[-1]['rel']:.2e}", flush=True)

    tau = 3
    _e2, ki2, ko2, p2 = one_structure(1500, tau, 2, 2, 0.0, 31)
    eps, lam, mu = 0.01, 0.8, 1.0
    S = states(tau)
    idx = {s: i for i, s in enumerate(S)}
    ident = []
    for dt in (0.008, 0.004, 0.002, 0.001):
        x = [0.0] * len(S)
        for a in range(tau + 1):
            x[idx[(a, tau - a, 0)]] = math.comb(tau, a) * (1 - eps) ** a * eps ** (tau - a)
        worst = 0.0
        for _n in range(int(20 / dt)):
            Phi = sum(x)
            worst = max(worst, abs(x[idx[(tau, 0, 0)]]
                                   - ((1 - eps) * p2.psi_tail(Phi)) ** tau))
            PhiA = sum(v for s, v in zip(S, x) if s[1] >= 1)
            h = lam * PhiA * p2.dpsi_tail(Phi) / p2.psi_tail(Phi)
            dx = [0.0] * len(S)
            for (a, b, c), v in zip(S, x):
                i = idx[(a, b, c)]
                if a:
                    dx[i] -= h * a * v; dx[idx[(a - 1, b + 1, c)]] += h * a * v
                if b:
                    dx[i] -= mu * b * v; dx[idx[(a, b - 1, c + 1)]] += mu * b * v
                if b >= 1:
                    dx[i] -= lam * v
            x = [q + dt * w for q, w in zip(x, dx)]
        ident.append({"dt": dt, "residual": worst})
        print(f"  dt={dt:.3f}: residual {worst:.3e}", flush=True)
    return {"bisection": bis, "identity": ident}


if __name__ == "__main__":
    data = {}
    print("finite size:"); data["finite_size"] = finite_size()
    print("threshold sweep:"); data["threshold"] = threshold_sweep()
    print("r_io scan:"); data["rio"] = rio_scan()
    print("trajectory:"); data["trajectory"] = trajectory()
    print("convergence:"); data["convergence"] = convergence_tables()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(data, f, indent=1)
    print("wrote", OUT)
