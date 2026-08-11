#!/usr/bin/env python3
"""Run the five verifications of section 3.1 and write every number to
data/verification.json. Nothing downstream recomputes anything: the figures and
the write-up both read this file, so a figure cannot disagree with the text.
"""
import json, math, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_closure import build, gillespie, Closure
from meanfield import MeanField
from audit_theory import Nmat, rho, C

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "data", "verification.json")

P = {(2, 2): .5, (3, 3): .5}
MS = (3, 4)
W = (1.0, 1.0)
TH = (1, 1)

def threshold_bp(P, ms, w, th):
    lo, hi = 1e-4, 6.0
    for _ in range(70):
        m = .5*(lo+hi)
        if rho(Nmat(P, ms, [m*x for x in w], th)) > 1: hi = m
        else: lo = m
    return .5*(lo+hi)

def jacobian_threshold(P, ms, w, th):
    """bisect on the largest real part of the closure's Jacobian at the
    disease-free state -- a route that shares no code with Nmat()."""
    def growth(lam):
        cl = Closure(P, ms, [lam*x for x in w], list(th)); cl.eps = 0.0
        n = cl.off[-1]
        x0 = [0.0]*(n+1)
        for a in range(cl.M):
            x0[cl.off[a]+cl.ix[a][(ms[a]-1, 0, 0)]] = 1.0
        f0 = cl.rhs(x0); h = 1e-7
        J = [[0.0]*(n) for _ in range(n)]
        for j in range(n):
            xp = list(x0); xp[j] += h
            fp = cl.rhs(xp)
            for q in range(n): J[q][j] = (fp[q]-f0[q])/h
        import numpy as _np
        return max(_np.linalg.eigvals(_np.array(J)).real)
    lo, hi = 1e-4, 6.0
    for _ in range(60):
        m = .5*(lo+hi)
        if growth(m) > 0: hi = m
        else: lo = m
    return .5*(lo+hi)

res = {}

# ---- 1. sum rule and the all-susceptible identity -------------------------
print("[1/5] sum rule and identity")
rows = []
for dt in (0.02, 0.01, 0.005, 0.0025):
    c = Closure(P, MS, [0.35, 0.35])
    idn, sm = c.invariants(0.01, 12.0, dt)
    rows.append({"dt": dt, "identity": idn, "sumrule": sm})
    print(f"    dt={dt}: identity {idn:.3e}  sum rule {sm:.3e}")
res["invariants"] = rows

# ---- 2. Jacobian of (2.5) versus rho(N) = 1 -------------------------------
print("[2/5] closure Jacobian versus rho(N)=1")
CASES = [
    ("uncorrelated", {(2, 2): .25, (2, 4): .25, (4, 2): .25, (4, 4): .25}, (3, 4), (1.0, 1.0), (1, 1)),
    ("positively correlated", {(2, 2): .5, (4, 4): .5}, (3, 4), (1.0, 1.0), (1, 1)),
    ("negatively correlated", {(2, 4): .5, (4, 2): .5}, (3, 4), (1.0, 1.0), (1, 1)),
    ("asymmetric rates and sizes", {(2, 2): .5, (4, 4): .5}, (3, 5), (1.0, 0.4), (1, 1)),
    ("pairwise limit", {(2, 3): .5, (4, 1): .5}, (2, 2), (1.0, 1.0), (1, 1)),
    ("three layers", {(2, 2, 2): .5, (3, 1, 2): .5}, (3, 4, 2), (1.0, .7, 1.3), (1, 1, 1)),
]
rows = []
for name, p, ms, w, th in CASES:
    a = jacobian_threshold(p, ms, w, th)
    b = threshold_bp(p, ms, w, th)
    rows.append({"case": name, "jacobian": a, "spectral": b, "rel": abs(a-b)/b})
    print(f"    {name:28s} {a:.8f}  {b:.8f}  rel {abs(a-b)/b:.1e}")
res["threshold_crosscheck"] = rows

# ---- 3. C: recursion versus single-group Monte Carlo ----------------------
print("[3/5] cascade C: recursion versus Monte Carlo")
def C_mc(m, lam, theta=1, trials=400000, seed=3):
    rng = random.Random(seed); tot = 0
    for _ in range(trials):
        i, s, new = 1, m-1, 0
        while i > 0 and s > 0:
            act = lam*s if i >= theta else 0.0
            T = act+i
            if T <= 0: break
            if rng.random() < act/T: i += 1; s -= 1; new += 1
            else: i -= 1
        tot += new
    return tot/trials
rows = []
for m, lam in ((2, 1.0), (3, 0.5), (3, 1.0), (4, 0.5), (4, 1.0), (5, 2.0), (6, 0.8)):
    a = C(m, lam); b = C_mc(m, lam)
    se = math.sqrt(max(b, 1e-9)/400000)*2
    rows.append({"m": m, "lam": lam, "recursion": a, "mc": b, "mc_se": se,
                 "naive": (m-1)*lam/(1+lam)})
    print(f"    m={m} lam={lam}: rec {a:.5f}  MC {b:.5f}+-{se:.5f}  (m-1)T {(m-1)*lam/(1+lam):.5f}")
res["cascade"] = rows

# ---- 4+5. closure and mean field against exact simulation ----------------
print("[4/5] trajectory, and [5/5] finite-size scaling")
lc = threshold_bp(P, MS, W, TH)
lam = 1.6*lc
betas = [lam*x for x in W]
eps, tmax, gdt = 0.02, 12.0, 0.25
res["setup"] = {"P": {str(k): v for k, v in P.items()}, "ms": list(MS), "w": list(W),
                "lambda_c": lc, "lambda": lam, "ratio": 1.6, "eps": eps}

cl_tr, _ = Closure(P, MS, betas).run(eps, tmax, 0.002)
mf_tr = MeanField(P, MS, betas).run(eps, tmax, 0.002)
grid = [round(j*gdt, 6) for j in range(int(tmax/gdt)+1)]
pick = lambda tr: {round(t, 6): (s, i) for t, s, i in tr}
clD, mfD = pick(cl_tr), pick(mf_tr)

def simulate(N0, runs, seed0):
    rng = random.Random(seed0)
    S, I = [], []
    for r in range(runs):
        N, layers = build(P, MS, N0, random.Random(seed0*31+r))
        tr = gillespie(N, layers, betas, tmax, gdt, eps, rng, list(TH))
        S.append([x[1] for x in tr][:len(grid)])
        I.append([x[2] for x in tr][:len(grid)])
    n = len(S)
    mS = [sum(p[j] for p in S)/n for j in range(len(grid))]
    mI = [sum(p[j] for p in I)/n for j in range(len(grid))]
    vS = [sum((p[j]-mS[j])**2 for p in S)/(n-1) for j in range(len(grid))]
    seS = [math.sqrt(v/n) for v in vS]
    vI = [sum((p[j]-mI[j])**2 for p in I)/(n-1) for j in range(len(grid))]
    seI = [math.sqrt(v/n) for v in vI]
    return mS, seS, mI, seI, n

mS, seS, mI, seI, n = simulate(6000, 400, 991)
res["trajectory"] = {"N": 6000, "runs": n, "t": grid,
                     "S_sim": mS, "S_sem": seS, "I_sim": mI, "I_sem": seI,
                     "S_closure": [clD[t][0] for t in grid],
                     "I_closure": [clD[t][1] for t in grid],
                     "S_mf": [mfD[t][0] for t in grid],
                     "I_mf": [mfD[t][1] for t in grid]}
print(f"    trajectory done (N=6000, {n} runs, lambda={lam:.4f}, lambda_c={lc:.4f})")

rows = []
for N0, runs in ((500, 600), (1000, 600), (2000, 800), (4000, 800), (8000, 500)):
    mS, seS, _, _, n = simulate(N0, runs, 1200+N0)
    dcl = sum(abs(mS[j]-clD[grid[j]][0]) for j in range(len(grid)))/len(grid)
    dmf = sum(abs(mS[j]-mfD[grid[j]][0]) for j in range(len(grid)))/len(grid)
    floor = sum(seS)/len(seS)*math.sqrt(2/math.pi)
    # RMS deviation with the sampling variance removed in quadrature
    def corrected(ref):
        m2 = sum((mS[j]-ref[grid[j]][0])**2 for j in range(len(grid)))/len(grid)
        v2 = sum(seS[j]**2 for j in range(len(grid)))/len(grid)
        return math.sqrt(max(m2-v2, 0.0))
    rows.append({"N": N0, "runs": n, "closure": dcl, "meanfield": dmf, "noise": floor,
                 "closure_corr": corrected(clD), "meanfield_corr": corrected(mfD)})
    print(f"    N={N0:5d} runs={n:4d}  closure {dcl:.5f} (corr {corrected(clD):.5f})"
          f"  mean field {dmf:.5f}  noise {floor:.5f}")
res["scaling"] = rows

os.makedirs(os.path.dirname(OUT), exist_ok=True)
json.dump(res, open(OUT, "w"), indent=1)
print("wrote", OUT)
