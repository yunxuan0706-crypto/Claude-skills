#!/usr/bin/env python3
"""The asymmetric-layer case was the one configuration whose closure deviation
sat above the noise floor at N=4000. Finite-size effect, or a real failure of
the closure for unequal rates and group sizes? Only the N-scaling tells them
apart: a finite-size effect keeps falling, a closure error plateaus.
"""
import math, random, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from verify_closure import build, gillespie, Closure
from audit_theory import Nmat, rho

P = {(2, 2): .5, (4, 4): .5}; ms = (3, 5); w = (1.0, 0.4); th = (1, 1)
lo, hi = 1e-3, 4.0
for _ in range(60):
    m = .5*(lo+hi)
    if rho(Nmat(P, ms, [m*x for x in w], th)) > 1: hi = m
    else: lo = m
lam = 1.6*.5*(lo+hi)
betas = [lam*x for x in w]
eps, tmax, gdt = 0.02, 12.0, 0.5
cl = Closure(P, ms, betas, list(th))
ode, _ = cl.run(eps, tmax, 0.005)
odeS = {round(t, 6): s for t, s in ode}
grid = [round(j*gdt, 6) for j in range(int(tmax/gdt)+1)]
print(f"asymmetric layers: m={ms}, w={w}, lambda={lam:.4f}")
print(f"{'N':>7} {'runs':>5} {'mean |dS|':>11} {'noise':>9} {'ratio':>7}")
for N0, runs in ((1000, 400), (2000, 400), (4000, 400), (8000, 250), (16000, 150)):
    rng = random.Random(31+N0); per = []
    for r in range(runs):
        N, layers = build(P, ms, N0, random.Random(777+N0*7+r))
        tr = gillespie(N, layers, betas, tmax, gdt, eps, rng, list(th))
        per.append([s for _, s in tr][:len(grid)])
    n = len(per)
    mean = [sum(p[j] for p in per)/n for j in range(len(grid))]
    var = [sum((p[j]-mean[j])**2 for p in per)/(n-1) for j in range(len(grid))]
    sem = [math.sqrt(v/n) for v in var]
    dev = sum(abs(mean[j]-odeS[grid[j]]) for j in range(len(grid)))/len(grid)
    fl = sum(sem)/len(sem)*math.sqrt(2/math.pi)
    print(f"{N0:>7} {n:>5} {dev:11.5f} {fl:9.5f} {dev/fl:7.2f}")
