#!/usr/bin/env python3
"""Audit the theory on the configurations the first checks never touched:
anti-correlated layer degrees, asymmetric rates and group sizes, and a layer
with theta = 2. Each case compares the full nonlinear closure with exact
simulation, and the closure's own nonlinear threshold with rho(N) = 1.
"""
import math, random, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from verify_closure import build, gillespie, Closure
from functools import lru_cache

def C(m, lam, theta=1):
    @lru_cache(None)
    def u(i, s):
        if s == 0 or i == 0: return 0.0
        act = lam*s if i >= theta else 0.0
        tot = act + i
        return (act/tot)*(1+u(i+1, s-1)) + (i/tot)*u(i-1, s)
    return u(1, m-1)

def Nmat(P, ms, betas, thetas):
    M = len(ms)
    kbar = [sum(p*k[a] for k, p in P.items()) for a in range(M)]
    A = [[0.0]*M for _ in range(M)]
    for a in range(M):
        for b in range(M):
            kab = sum(p*k[a]*k[b] for k, p in P.items())
            A[a][b] = C(ms[b], betas[b], thetas[b])*(kab/kbar[a] - (1 if a == b else 0))
    return A

def rho(A):
    n = len(A)
    v = [1.0]*n
    for _ in range(3000):
        w = [sum(A[i][j]*v[j] for j in range(n)) for i in range(n)]
        nv = math.sqrt(sum(x*x for x in w))
        if nv == 0: return 0.0
        v = [x/nv for x in w]
    return nv

def ode_threshold(P, ms, thetas, w, eps=1e-9, tmax=400.0, dt=0.02):
    """bisect on lambda using the nonlinear closure's own final size"""
    def macro(L):
        cl = Closure(P, ms, [L*wa for wa in w], thetas)
        _, y = cl.run(eps, tmax, dt)
        return y[-1] > 1e-3            # R(inf) macroscopic?
    lo, hi = 1e-3, 4.0
    if not macro(hi): return float('nan')
    for _ in range(45):
        mid = .5*(lo+hi)
        if macro(mid): hi = mid
        else: lo = mid
    return .5*(lo+hi)

def bp_threshold(P, ms, thetas, w):
    lo, hi = 1e-3, 4.0
    for _ in range(60):
        mid = .5*(lo+hi)
        if rho(Nmat(P, ms, [mid*wa for wa in w], thetas)) > 1: hi = mid
        else: lo = mid
    return .5*(lo+hi)

if __name__ == "__main__":
    CASES = [
        ("anti-correlated degrees", {(2, 4): .5, (4, 2): .5}, (3, 4), (1.0, 1.0), (1, 1)),
        ("asymmetric rates+sizes",  {(2, 2): .5, (4, 4): .5}, (3, 5), (1.0, 0.4), (1, 1)),
        ("theta = (1,2)",           {(2, 2): .5, (3, 3): .5}, (3, 4), (1.0, 1.0), (1, 2)),
        ("three layers",            {(2, 2, 2): .5, (3, 1, 2): .5}, (3, 4, 2), (1.0, .7, 1.3), (1, 1, 1)),
    ]

    print("A. nonlinear closure threshold  vs  rho(N)=1")
    print(f"{'case':26s} {'ODE final-size':>15s} {'rho(N)=1':>10s} {'rel':>9s}")
    for name, P, ms, w, th in CASES:
        a = ode_threshold(P, ms, th, w); b = bp_threshold(P, ms, th, w)
        rel = abs(a-b)/b if b > 0 and a == a else float('nan')
        print(f"{name:26s} {a:15.5f} {b:10.5f} {rel:9.1e}")

    print("\nB. full nonlinear closure vs exact simulation (N = 4000, 300 runs)")
    print(f"{'case':26s} {'mean |dS|':>10s} {'noise':>9s} {'verdict':>10s}")
    for name, P, ms, w, th in CASES:
        lam = 1.6*bp_threshold(P, ms, th, w)
        if not (lam == lam) or lam > 3: 
            print(f"{name:26s} {'--- no finite threshold ---':>32s}"); continue
        betas = [lam*wa for wa in w]
        eps, tmax, gdt, N0, runs = 0.02, 12.0, 0.5, 4000, 300
        cl = Closure(P, ms, betas, th)
        ode, _ = cl.run(eps, tmax, 0.005)
        odeS = {round(t, 6): s for t, s, _ in ode}
        grid = [round(j*gdt, 6) for j in range(int(tmax/gdt)+1)]
        rng = random.Random(99)
        per = []
        for r in range(runs):
            N, layers = build(P, ms, N0, random.Random(5150+r))
            tr = gillespie(N, layers, betas, tmax, gdt, eps, rng, list(th))
            per.append([x[1] for x in tr][:len(grid)])
        n = len(per)
        mean = [sum(p[j] for p in per)/n for j in range(len(grid))]
        var = [sum((p[j]-mean[j])**2 for p in per)/(n-1) for j in range(len(grid))]
        sem = [math.sqrt(v/n) for v in var]
        dev = sum(abs(mean[j]-odeS[grid[j]]) for j in range(len(grid)))/len(grid)
        floor = sum(sem)/len(sem)*math.sqrt(2/math.pi)
        print(f"{name:26s} {dev:10.5f} {floor:9.5f} {'OK' if dev < 3*floor else 'CHECK':>10s}")
