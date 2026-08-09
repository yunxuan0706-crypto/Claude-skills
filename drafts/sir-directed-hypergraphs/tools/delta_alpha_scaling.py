"""Does the breaking signal scale with Delta-alpha at all?

Result (N=200, K=8 structures, 6000 runs/point, lambda = 1.3 lambda_c):

    Delta-alpha   0.000   0.036   0.075   0.108
    <dPi>        -0.0025 -0.0021 -0.0022 -0.0017
    sigma          0.97    0.96    0.90    0.58

Flat. The Delta-alpha = 0 baseline is indistinguishable from Delta-alpha =
0.108, so there is no dependence at all; the small common offset is a
pairing systematic, not a structural effect.

In the SAME protocol, order asymmetry (tau=3, eta=1) gives
<dPi> = +0.0384 +/- 0.0021, i.e. 18 sigma. The protocol detects breaking
when breaking is there.

Conclusion: direction-symmetry breaking in this model is driven by the
reversal-odd part of the BRANCHING structure (order asymmetry, bi-degree
sequence asymmetry), not by the overlap polarity. Being reversal-odd is
necessary but nowhere near sufficient.
"""
import math, random, statistics as st, sys
sys.path.insert(0,'tools')
from delta_alpha_breaking import mirror_structure, Overlap, drive, sizes, summarise
from feasibility_check import index, gillespie_fast, joint_of, mf_threshold

Nh, Mh, tau, eta, K, RUNS = 100, 75, 2, 2, 8, 6000
N = 2*Nh
probe = mirror_structure(Nh, Mh, tau, eta, random.Random(1))
jc,_,_ = joint_of(N, probe); lc = mf_threshold(jc, tau)
lam = 1.3*lc
print(f"N={N} M={2*Mh} K={K} runs={RUNS} lam=1.3*lc={lam:.4f}\n")

print("--- Delta-alpha scaling (each row = pooled over K structures) ---")
print(f"{'proposals':>10} {'mean Delta':>11} {'<dPi>':>11} {'+/-':>9} {'sigma':>7}")
for props in (0, 8000, 30000, 90000):
    ds, dpis = [], []
    for k in range(K):
        edges = mirror_structure(Nh, Mh, tau, eta, random.Random(500+13*k))
        ov = Overlap(edges, tau, eta)
        if props: drive(edges, ov, random.Random(900+7*k), props, None)
        ds.append(ov.delta())
        edgesR = [(H,T) for T,H in edges]
        SA = sizes(N, edges , lam, RUNS, 200000+900*k)
        SB = sizes(N, edgesR, lam, RUNS, 200000+900*k)
        pA,_ = summarise(SA, 0.02); pB,_ = summarise(SB, 0.02)
        dpis.append(pA-pB)
    md = st.mean(dpis); sem = math.sqrt(st.variance(dpis)/K)
    print(f"{props:10d} {st.mean(ds):+11.4f} {md:+11.5f} {sem:9.5f} {abs(md)/sem:7.2f}")

print("\n--- positive control in the SAME protocol: order asymmetry tau=3,eta=1 ---")
dpis=[]
for k in range(K):
    rng=random.Random(3000+11*k); edges=[]
    for _ in range(2*Mh):
        nd=rng.sample(range(N),4); edges.append((set(nd[:3]),set(nd[3:])))
    edgesR=[(H,T) for T,H in edges]
    jc2,_,_=joint_of(N,edges); lc2=mf_threshold(jc2,3); l2=1.3*lc2
    SA=sizes(N,edges ,l2,RUNS,400000+900*k)
    SB=sizes(N,edgesR,l2,RUNS,400000+900*k)
    pA,_=summarise(SA,0.02); pB,_=summarise(SB,0.02)
    dpis.append(pA-pB)
md=st.mean(dpis); sem=math.sqrt(st.variance(dpis)/K)
print(f"  <dPi> = {md:+.5f} +/- {sem:.5f}  ({abs(md)/sem:.2f} sigma)  -> protocol CAN detect breaking")
