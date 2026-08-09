"""Does the EBCM's overshoot at lambda_c actually shrink with N?

Section III.E attributes the +18% overshoot at threshold to the closure being
an N->infinity theory. That is a testable claim, not an explanation: if the
overshoot were a closure error it would persist at large N.
"""
import math, os, random, statistics as st, sys
sys.path.insert(0, "/home/user/Claude-skills/drafts/sir-directed-hypergraphs/tools")
from ebcm_directed import PGF, bidegree, config_model, degrees, ebcm, gillespie_traj
from feasibility_check import index

eps, reps = 0.005, 3
print(f"eps={eps}, {reps} structures per N, at lambda = lambda_c exactly")
print("   N   runs |  sim R(inf)        EBCM      overshoot")
for N, runs in ((1000, 900), (2000, 700), (4000, 500), (8000, 400), (16000, 260)):
    over = []
    for r in range(reps):
        rng = random.Random(9100 + 37 * r + N)
        kin, kout = bidegree(N, 2, 2, 2, 0.0, rng)
        edges = config_model(kin, kout, 2, 2, rng)
        ki, ko = degrees(N, edges)
        p = PGF(ki, ko)
        lam = 1.0 / (2 * p.kappa - 1)
        tails, heads, tails_of = index(N, edges)
        vals = []
        for s in range(runs):
            rr = random.Random(4_100_000 + 7919 * s + 13 * r + N)
            sd = [v for v in range(N) if rr.random() < eps]
            vals.append(gillespie_traj(N, tails, heads, tails_of, lam, 1, rr,
                                       sd, [400.0])[-1][2])
        m = st.mean(vals)
        e = ebcm(p, 2, lam, eps=eps, tmax=250.0, dt=0.01)[-1][3]
        over.append((e - m) / m)
        if r == 0:
            first = (m, math.sqrt(st.variance(vals) / runs), e)
    mo = st.mean(over)
    so = math.sqrt(st.variance(over) / len(over)) if len(over) > 1 else 0.0
    print(f"{N:6d} {runs:5d} | {first[0]:.5f}+/-{first[1]:.5f}  {first[2]:.5f}"
          f"   {100*mo:+6.1f}% +/- {100*so:.1f}%")
