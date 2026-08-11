#!/usr/bin/env python3
"""Recompute the scaling block with a bootstrap interval on the deviation.

The deviation statistic is mean_t |mean_sim(t) - theory(t)| over a time grid.
The pointwise sampling error is easy (the standard error of mean_sim at each t),
but the deviations are strongly autocorrelated in time (lag-1 ~ 0.95), so the
grid carries far fewer than 49 independent pieces of information and a pointwise
"noise floor" says nothing about the uncertainty of the average itself.
Resampling whole runs with replacement handles the time correlation for free and
gives an interval on the statistic that is actually being compared.
"""
import json, math, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_closure import build, gillespie, Closure
from meanfield import MeanField
from audit_theory import Nmat, rho

HERE = os.path.dirname(os.path.abspath(__file__))
F = os.path.join(HERE, "..", "data", "verification.json")
d = json.load(open(F))
s = d["setup"]
P = {(2, 2): .5, (3, 3): .5}
MS, W, TH = tuple(s["ms"]), tuple(s["w"]), (1, 1)
betas = [s["lambda"]*x for x in W]
eps, tmax, gdt = s["eps"], 12.0, 0.25

cl, _ = Closure(P, MS, betas).run(eps, tmax, 0.002)
mf = MeanField(P, MS, betas).run(eps, tmax, 0.002)
grid = [round(j*gdt, 6) for j in range(int(tmax/gdt)+1)]
clD = {round(t, 6): S for t, S, _ in cl}
mfD = {round(t, 6): S for t, S, _ in mf}
ref_cl = [clD[t] for t in grid]
ref_mf = [mfD[t] for t in grid]

def boot(per, ref, B=2000, seed=17):
    """bootstrap over runs -> percentile interval on mean_t |mean_sim - ref|"""
    rng = random.Random(seed); n = len(per); T = len(ref)
    stats = []
    for _ in range(B):
        idx = [rng.randrange(n) for _ in range(n)]
        m = [sum(per[i][j] for i in idx)/n for j in range(T)]
        stats.append(sum(abs(m[j]-ref[j]) for j in range(T))/T)
    stats.sort()
    return stats[int(.025*B)], stats[int(.975*B)]

rows = []
print(f"{'N':>7} {'runs':>5} {'closure [95% CI]':>28} {'mean field [95% CI]':>28}")
for N0, runs in ((500, 600), (1000, 600), (2000, 800), (4000, 800), (8000, 500)):
    rng = random.Random(1200+N0)
    per = []
    for r in range(runs):
        N, layers = build(P, MS, N0, random.Random((1200+N0)*31+r))
        tr = gillespie(N, layers, betas, tmax, gdt, eps, rng, list(TH))
        per.append([x[1] for x in tr][:len(grid)])
    n = len(per); T = len(grid)
    mS = [sum(p[j] for p in per)/n for j in range(T)]
    seS = [math.sqrt(sum((p[j]-mS[j])**2 for p in per)/(n-1)/n) for j in range(T)]
    dcl = sum(abs(mS[j]-ref_cl[j]) for j in range(T))/T
    dmf = sum(abs(mS[j]-ref_mf[j]) for j in range(T))/T
    lcl, ucl = boot(per, ref_cl)
    lmf, umf = boot(per, ref_mf)
    floor = sum(seS)/T*math.sqrt(2/math.pi)
    rows.append({"N": N0, "runs": n, "closure": dcl, "closure_lo": lcl, "closure_hi": ucl,
                 "meanfield": dmf, "meanfield_lo": lmf, "meanfield_hi": umf,
                 "noise": floor})
    print(f"{N0:>7} {n:>5}   {dcl:.5f} [{lcl:.5f},{ucl:.5f}]   {dmf:.5f} [{lmf:.5f},{umf:.5f}]")

d["scaling"] = rows
json.dump(d, open(F, "w"), indent=1)
print("updated", F)
