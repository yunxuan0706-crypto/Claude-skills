#!/usr/bin/env python3
"""Finite-size scaling of the closure error, with the Monte-Carlo noise floor.

The first run showed the error rising again at the largest N. Without an error
bar that is unreadable: it could be a real breakdown of the closure or simply
the run count. Report both the deviation and the sampling error on it.
"""
import random, math, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from verify_closure import build, gillespie, Closure

P = {(2, 2): .5, (3, 3): .5}; ms = (3, 4); lam = 0.35
betas = [lam, lam]; eps, tmax, gdt = 0.02, 12.0, 0.25
cl = Closure(P, ms, betas)
ode, _ = cl.run(eps, tmax, 0.005)
odeS = {round(t, 6): s for t, s in ode}
grid = [round(j*gdt, 6) for j in range(int(tmax/gdt)+1)]
ref = [odeS[t] for t in grid]

print(f"lambda={lam}, m={ms}, eps={eps};  mean |dS| over t<=12 with its sampling error")
print(f"{'N':>7} {'runs':>5} {'mean |dS|':>11} {'noise floor':>12} {'above noise?':>13}")
for N0, runs in ((500, 500), (1000, 500), (2000, 500), (4000, 500), (8000, 300)):
    rng = random.Random(11 + N0)
    per = []
    for r in range(runs):
        N, layers = build(P, ms, N0, random.Random(4242 + N0*17 + r))
        tr = gillespie(N, layers, betas, tmax, gdt, eps, rng)
        per.append([s for _, s in tr][:len(grid)])
    n = len(per)
    mean = [sum(p[j] for p in per)/n for j in range(len(grid))]
    var = [sum((p[j]-mean[j])**2 for p in per)/(n-1) for j in range(len(grid))]
    sem = [math.sqrt(v/n) for v in var]
    dev = sum(abs(mean[j]-ref[j]) for j in range(len(grid)))/len(grid)
    floor = sum(sem)/len(sem)*math.sqrt(2/math.pi)     # E|noise| for a centred normal
    print(f"{N0:>7} {n:>5} {dev:11.5f} {floor:12.5f} {'yes' if dev > 2*floor else 'NO':>13}")
