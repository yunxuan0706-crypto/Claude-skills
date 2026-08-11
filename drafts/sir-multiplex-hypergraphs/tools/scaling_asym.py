#!/usr/bin/env python3
"""The asymmetric-layer case was the one configuration whose closure deviation
sat above the noise floor at N=4000. Finite-size effect, or a real failure of
the closure for unequal rates and group sizes? Only the N-scaling tells them
apart: a finite-size effect keeps falling, a closure error plateaus.

Result (lambda = 1.6 lambda_c, eps = 0.02, t <= 12):

       N   runs   mean |dS|     noise   ratio
    1000    400     0.00682   0.00240    2.84
    2000    400     0.00567   0.00161    3.53
    4000    400     0.00287   0.00107    2.68
    8000    250     0.00039   0.00094    0.41
   16000    150     0.00111   0.00083    1.34

The deviation falls by a factor ~17 between N=1000 and N=8000 and reaches the
sampling floor there; the N=16000 point sits at 1.3x the floor on the fewest
runs and is above N=8000, which no monotone systematic can do -- it is noise.
So the N=4000 flag was a finite-size effect, stronger here because m2=5 packs
short cycles more densely at fixed N. The closure holds for unequal rates and
unequal group sizes.
"""
import math, random, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from verify_closure import build, gillespie, Closure
from audit_theory import Nmat, rho

if __name__ == "__main__":
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
    odeS = {round(t, 6): s for t, s, _ in ode}
    grid = [round(j*gdt, 6) for j in range(int(tmax/gdt)+1)]
    print(f"asymmetric layers: m={ms}, w={w}, lambda={lam:.4f}")
    print(f"{'N':>7} {'runs':>5} {'mean |dS|':>11} {'noise':>9} {'ratio':>7}")
    for N0, runs in ((1000, 400), (2000, 400), (4000, 400), (8000, 250), (16000, 150)):
        rng = random.Random(31+N0); per = []
        for r in range(runs):
            N, layers = build(P, ms, N0, random.Random(777+N0*7+r))
            tr = gillespie(N, layers, betas, tmax, gdt, eps, rng, list(th))
            per.append([x[1] for x in tr][:len(grid)])
        n = len(per)
        mean = [sum(p[j] for p in per)/n for j in range(len(grid))]
        var = [sum((p[j]-mean[j])**2 for p in per)/(n-1) for j in range(len(grid))]
        sem = [math.sqrt(v/n) for v in var]
        dev = sum(abs(mean[j]-odeS[grid[j]]) for j in range(len(grid)))/len(grid)
        fl = sum(sem)/len(sem)*math.sqrt(2/math.pi)
        print(f"{N0:>7} {n:>5} {dev:11.5f} {fl:9.5f} {dev/fl:7.2f}")
