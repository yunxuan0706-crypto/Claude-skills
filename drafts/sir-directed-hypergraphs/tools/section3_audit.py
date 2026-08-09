#!/usr/bin/env python3
"""Independent audit of every Section III number.

Re-running the same code proves nothing. Each check here reaches the quantity
by a route that does not share the derivation being tested.

 A  Final size for tau=1 in closed form. The ODEs give dPhi/dt = -beta*x_I and
    dx_R/dt = mu*x_I, so x_R = (1/lambda)(1-Phi) exactly, and at t=inf with
    x_I=0,
        Phi_inf = [1 + lambda*(1-eps)*psi_tail(Phi_inf)] / (1 + lambda).
    Solving that scalar fixed point uses no integrator. Compare with R(inf)
    from the full RK4 run.

 B  R_0 by direct Monte Carlo. Draw a node the way a newly infected node is
    reached -- through a head slot, hence k_in-biased -- give it an Exp(mu)
    infectious period, and race an Exp(beta) clock on every (out-hyperedge,
    head) pair. The mean offspring count must equal tau*kappa*beta/(beta+mu).
    This tests the branching reading of eq (35) and the identity (36) at once,
    touching neither the ODE nor the final-size algebra.

 C  Mean-field threshold by bisection on its own integrated final size, the
    same way the EBCM threshold was checked. lambda_c^MF = 1/(tau*kappa) was
    until now only asserted.

 D  Generator hygiene. config_model repairs degenerate edges by swapping slots
    and then DROPS whatever is still degenerate, so the realised M and the
    realised bi-degree sequence can differ from the requested ones. Measure the
    loss; the PGFs are built from the realised sequence, but a large loss would
    distort the r_io targeting.

 E  Seeding mismatch. The ODE starts from exactly eps; the simulation seeds
    each node with probability eps, so the seed count fluctuates by ~sqrt(eps*N).
    Compare Bernoulli seeding against a fixed seed count.

 F  Integrator convergence of the reported trajectory (dt -> dt/2).

 G  Robustness of the fitted exponent and of the mean-field plateau to which
    points enter the fit.

 H  The fast Gillespie sampler against the naive reference implementation.
"""
import math
import os
import random
import statistics as st
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ebcm_directed import (PGF, bidegree, config_model, degrees, ebcm,
                           gillespie_traj, meanfield)
from feasibility_check import index, gillespie_fast, gillespie_ref


def rule(t):
    print("\n" + t + "\n" + "-" * len(t))


# ------------------------------------------------------------------- A ------
def check_A():
    rule("A  tau=1 final size: closed-form fixed point vs integrated ODE")
    for seed, md, lam_r in ((11, 3, 1.6), (12, 2, 2.5), (13, 4, 1.3)):
        rng = random.Random(seed)
        kin, kout = bidegree(4000, 1, 1, md, 0.0, rng)
        edges = config_model(kin, kout, 1, 1, rng)
        ki, ko = degrees(4000, edges)
        p = PGF(ki, ko)
        lam = lam_r / (p.kappa - 1)
        eps = 0.01
        # scalar fixed point, no integrator
        Phi = 1.0
        for _ in range(20000):
            nxt = (1 + lam * (1 - eps) * p.psi_tail(Phi)) / (1 + lam)
            if abs(nxt - Phi) < 1e-15:
                Phi = nxt
                break
            Phi = nxt
        R_closed = 1 - (1 - eps) * p.psi_in(Phi)
        R_ode = ebcm(p, 1, lam, eps=eps, tmax=400.0, dt=0.002)[-1][3]
        print(f"  kappa={p.kappa:.4f} lam={lam:.4f}: closed {R_closed:.9f}"
              f"  ODE {R_ode:.9f}  |diff| {abs(R_closed-R_ode):.2e}")


# ------------------------------------------------------------------- B ------
def check_B(N=4000, tau=2, eta=2, md=2, trials=400000, seed=77):
    rule("B  R_0 by direct Monte Carlo vs tau*kappa*beta/(beta+mu)")
    rng = random.Random(seed)
    kin, kout = bidegree(N, tau, eta, md, 0.3, rng)
    edges = config_model(kin, kout, tau, eta, rng)
    ki, ko = degrees(N, edges)
    p = PGF(ki, ko)
    # head slots: reaching a node the way an infection reaches it
    head_slots = [v for v in range(N) for _ in range(ki[v])]
    out_of = [[] for _ in range(N)]
    for i, (T, H) in enumerate(edges):
        for v in T:
            out_of[v].append(i)
    for lam in (0.3, 0.8, 2.0):
        beta, mu = lam, 1.0
        tot = 0
        for _ in range(trials):
            v = head_slots[rng.randrange(len(head_slots))]
            life = rng.expovariate(mu)
            for f in out_of[v]:
                for _w in edges[f][1]:
                    if rng.expovariate(beta) < life:
                        tot += 1
        meas = tot / trials
        pred = tau * p.kappa * beta / (beta + mu)
        se = math.sqrt(meas / trials) * 2      # crude; offspring is over-dispersed
        print(f"  lam={lam:.2f}: measured R_0 = {meas:.4f}   predicted {pred:.4f}"
              f"   rel {abs(meas-pred)/pred:.2e}")
    # the two counting forms of eq (36)
    a = tau * p.kappa
    b = eta * sum(x * y for x, y in zip(ki, ko)) / sum(ki)
    print(f"  eq (36): tau<kk>/<ko> = {a:.10f}   eta<kk>/<ki> = {b:.10f}"
          f"   rel {abs(a-b)/a:.1e}")


# ------------------------------------------------------------------- C ------
def check_C(N=1500, tau=2, eta=2, md=2, seed=23):
    rule("C  mean-field threshold by bisection vs 1/(tau*kappa)")
    rng = random.Random(seed)
    kin, kout = bidegree(N, tau, eta, md, 0.0, rng)
    edges = config_model(kin, kout, tau, eta, rng)
    ki, ko = degrees(N, edges)
    p = PGF(ki, ko)
    pred = 1.0 / (tau * p.kappa)
    for tmax in (400.0, 1600.0, 6400.0):
        lo, hi = 1e-3, 50.0
        while hi - lo > 1e-5 * hi:
            mid = math.sqrt(lo * hi)
            sol = meanfield(N, edges, ki, ko, tau, mid, eps=1e-12,
                            tmax=tmax, dt=0.05)
            if sol[-1][3] > 1e-8:
                hi = mid
            else:
                lo = mid
        meas = math.sqrt(lo * hi)
        print(f"  tmax={tmax:6.0f}: measured {meas:.6f}  predicted {pred:.6f}"
              f"  rel {abs(meas-pred)/pred:.2e}")


# ------------------------------------------------------------------- D ------
def check_D():
    rule("D  generator hygiene: how much does edge repair actually discard?")
    for N, tau, eta, md, corr in ((500, 2, 2, 2, 0.0), (3000, 2, 2, 2, 0.0),
                                  (3000, 2, 2, 2, 0.9), (3000, 2, 2, 2, -0.9),
                                  (4000, 3, 1, 2, 0.0)):
        rng = random.Random(4242)
        kin, kout = bidegree(N, tau, eta, md, corr, rng)
        M_req = sum(kout) // tau
        edges = config_model(kin, kout, tau, eta, rng)
        ki, ko = degrees(N, edges)
        lost = 1 - len(edges) / M_req
        dko = sum(abs(a - b) for a, b in zip(ko, kout)) / max(1, sum(kout))
        dki = sum(abs(a - b) for a, b in zip(ki, kin)) / max(1, sum(kin))
        print(f"  N={N} tau={tau} eta={eta} corr={corr:+.1f}: "
              f"edges kept {len(edges)}/{M_req} (lost {lost:.2%}), "
              f"|dk_out|/sum {dko:.2%}, |dk_in|/sum {dki:.2%}")


# ------------------------------------------------------------------- E ------
def check_E(N=6000, tau=2, eta=2, md=2, runs=400, eps=0.01, seed=4000):
    rule("E  Bernoulli seeding vs a fixed seed count (the ODE assumes exactly eps)")
    rng = random.Random(seed)
    kin, kout = bidegree(N, tau, eta, md, 0.0, rng)
    edges = config_model(kin, kout, tau, eta, rng)
    ki, ko = degrees(N, edges)
    p = PGF(ki, ko)
    lam = 1.6 / (tau * p.kappa - 1)
    tails, heads, tails_of = index(N, edges)
    k = int(round(eps * N))
    for mode in ("bernoulli", "fixed"):
        vals = []
        for s in range(runs):
            rr = random.Random(880000 + 13 * s)
            if mode == "bernoulli":
                sd = [v for v in range(N) if rr.random() < eps]
            else:
                sd = rr.sample(range(N), k)
            vals.append(gillespie_traj(N, tails, heads, tails_of, lam, 1, rr,
                                       sd, [400.0])[-1][2])
        m = st.mean(vals)
        se = math.sqrt(st.variance(vals) / runs)
        print(f"  {mode:9s}: R(inf) = {m:.5f} +/- {se:.5f}")


# ------------------------------------------------------------------- F ------
def check_F(N=6000, tau=2, eta=2, md=2, seed=4000):
    rule("F  integrator convergence of the reported trajectory")
    rng = random.Random(seed)
    kin, kout = bidegree(N, tau, eta, md, 0.0, rng)
    edges = config_model(kin, kout, tau, eta, rng)
    ki, ko = degrees(N, edges)
    p = PGF(ki, ko)
    lam = 1.6 / (tau * p.kappa - 1)
    prev = None
    for dt in (0.02, 0.01, 0.005, 0.0025):
        sol = ebcm(p, tau, lam, eps=0.01, tmax=14.0, dt=dt)
        R = sol[-1][3]
        d = "" if prev is None else f"   |change| {abs(R-prev):.2e}"
        print(f"  dt={dt:7.4f}: R(14) = {R:.10f}{d}")
        prev = R


# ------------------------------------------------------------------- G ------
def check_G():
    rule("G  are the fitted exponent and the mean-field plateau robust?")
    import json
    d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "data", "section3.json")))
    fs = d["finite_size"]
    N = [r["N"] for r in fs]

    def fit_pow(idx):
        lx = [math.log(N[i]) for i in idx]
        ly = [math.log(fs[i]["ebcm"]) for i in idx]
        w = [(fs[i]["ebcm"] / fs[i]["ebcm_sem"]) ** 2 for i in idx]
        Sw = sum(w); Sx = sum(a * b for a, b in zip(w, lx))
        Sy = sum(a * b for a, b in zip(w, ly))
        Sxx = sum(a * b * b for a, b in zip(w, lx))
        Sxy = sum(a * b * c for a, b, c in zip(w, lx, ly))
        det = Sw * Sxx - Sx * Sx
        return (Sw * Sxy - Sx * Sy) / det, math.sqrt(Sw / det)

    for label, idx in (("all four", [0, 1, 2, 3]), ("drop N=500", [1, 2, 3]),
                       ("drop N=4000", [0, 1, 2]), ("two largest", [2, 3])):
        s, e = fit_pow(idx)
        print(f"  EBCM exponent [{label:12s}] = {s:+.3f} +/- {e:.3f}")

    def plateau(idx):
        xs = [1 / N[i] for i in idx]
        ys = [fs[i]["mf"] for i in idx]
        n = len(xs); sx = sum(xs); sy = sum(ys)
        sxx = sum(x * x for x in xs); sxy = sum(x * y for x, y in zip(xs, ys))
        b = (n * sxy - sx * sy) / (n * sxx - sx * sx)
        return (sy - b * sx) / n

    for label, idx in (("all four", [0, 1, 2, 3]), ("drop N=500", [1, 2, 3]),
                       ("three largest", [1, 2, 3]), ("two largest", [2, 3])):
        print(f"  MF plateau  [{label:12s}] = {plateau(idx):.4f}")


# ------------------------------------------------------------------- H ------
def check_H(N=120, tau=2, eta=2, M=90, runs=3000, lam=1.2, seed=5):
    rule("H  fast sampler vs the naive reference implementation")
    rng = random.Random(seed)
    edges = []
    while len(edges) < M:
        nd = rng.sample(range(N), tau + eta)
        edges.append((set(nd[:tau]), set(nd[tau:])))
    tails, heads, tails_of = index(N, edges)
    a, b = [], []
    for s in range(runs):
        sd = [random.Random(9000 + s).randrange(N)]
        a.append(gillespie_fast(N, tails, heads, tails_of, lam, 1,
                                random.Random(31 * s + 1), sd))
        b.append(gillespie_ref(edges, N, lam, 1, random.Random(31 * s + 2), sd))
    ma, mb = st.mean(a), st.mean(b)
    se = math.sqrt(st.variance(a) / runs + st.variance(b) / runs)
    print(f"  fast {ma:.5f}   reference {mb:.5f}   diff {ma-mb:+.5f}"
          f"  ({abs(ma-mb)/se:.2f} sigma over {runs} runs)")


if __name__ == "__main__":
    check_A(); check_C(); check_D(); check_F(); check_G(); check_H()
    check_E(); check_B()
