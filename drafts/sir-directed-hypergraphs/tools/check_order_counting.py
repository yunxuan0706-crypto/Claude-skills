#!/usr/bin/env python3
"""Which states of the concurrent EBCM are first order near the disease-free state?

Section III.C linearises around x_(tau,0,0)=1. The text originally said states
with b>=2 are second order, which misclassifies x_(tau-2,1,1): it carries b=1
but needs two tail members to have been infected. The correct grading is by
b+c, the number of tail members EVER infected.

Measured by integrating from two seed sizes and reading off the power of eps:
a ratio of 10 between eps=1e-4 and 1e-5 means first order, 100 second, and so
on. Result: the order equals b+c exactly, with (tau-2,1,1) coming out second
order despite b=1. The linearisation itself is unaffected -- Phi_A's
first-order part is p either way -- but the stated criterion was wrong.
"""
import math, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ebcm_directed import PGF, bidegree, config_model, degrees, states

def main(tau=3, N=1200, ratio=0.9, T=3.0, dt=0.001, seed=7):
    rng = random.Random(seed)
    kin, kout = bidegree(N, tau, 2, 2, 0.0, rng)
    E = config_model(kin, kout, tau, 2, rng)
    ki, ko = degrees(N, E)
    p = PGF(ki, ko)
    lam = ratio / (tau * p.kappa - 1)
    beta, mu = lam, 1.0
    S = states(tau)
    idx = {s: i for i, s in enumerate(S)}

    def run(eps):
        x = [0.0] * len(S)
        for a in range(tau + 1):
            x[idx[(a, tau - a, 0)]] = math.comb(tau, a) * (1 - eps) ** a * eps ** (tau - a)
        for _ in range(int(T / dt)):
            Phi = sum(x)
            PhiA = sum(v for s, v in zip(S, x) if s[1] >= 1)
            h = beta * PhiA * p.dpsi_tail(Phi) / p.psi_tail(Phi)
            dx = [0.0] * len(S)
            for (a, b, c), v in zip(S, x):
                i = idx[(a, b, c)]
                if a:
                    dx[i] -= h * a * v; dx[idx[(a - 1, b + 1, c)]] += h * a * v
                if b:
                    dx[i] -= mu * b * v; dx[idx[(a, b - 1, c + 1)]] += mu * b * v
                if b >= 1:
                    dx[i] -= beta * v
            x = [q + dt * w for q, w in zip(x, dx)]
        return x

    x1, x2 = run(1e-4), run(1e-5)
    print(f"tau={tau}, lambda={lam:.4f} ({ratio}*lambda_c)")
    print(" state     b+c   x(1e-4)     x(1e-5)     ratio    order")
    bad = 0
    for s in S:
        a, b, c = s
        if (a, b, c) == (tau, 0, 0):
            continue
        v1, v2 = x1[idx[s]], x2[idx[s]]
        if v1 <= 0 or v2 <= 0:
            continue
        order = math.log10(v1 / v2)
        ok = abs(order - (b + c)) < 0.05
        bad += not ok
        print(f" ({a},{b},{c})      {b+c}   {v1:.3e}  {v2:.3e}  {v1/v2:8.1f}   "
              f"{order:.2f}{'' if ok else '   <-- != b+c'}")
    print(f"\nstates whose measured order differs from b+c: {bad}")

if __name__ == "__main__":
    main()
