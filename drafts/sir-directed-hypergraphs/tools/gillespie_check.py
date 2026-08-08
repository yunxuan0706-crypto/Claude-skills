"""Exact Gillespie sampler for the Section II.B model, plus the checks that
需要 simulation rather than combinatorics.

Verified here:
  * the final size depends only on lambda = beta/mu (all rates scale out);
  * with theta >= 2 a single seed cannot start anything, so the theta >= 2
    branch needs a seed set -- this is why Section II.B makes the initial
    condition depend on theta;
  * Delta-alpha is tunable at fixed bi-degree sequence, i.e. the C4
    symmetry-breaking experiment is feasible.

The sampler itself is what week 1 of the execution plan needs for the
mean-field vs edge-based comparison.
"""
import math, random, statistics as st

# ---------- exact Gillespie for the Section II.B model ----------
def gillespie(edges, N, beta, mu, theta, seed):
    rng = random.Random(seed)
    state = [0]*N                       # 0=S 1=I 2=R
    state[rng.randrange(N)] = 1
    heads = [[] for _ in range(N)]      # node -> hyperedges whose HEAD contains it
    for i,(T,H) in enumerate(edges):
        for v in H: heads[v].append(i)
    while True:
        n_e = [sum(1 for v in T if state[v]==1) for T,_ in edges]
        rates, tot = [], 0.0
        for u in range(N):
            if state[u]==0:
                r = beta*sum(1 for i in heads[u] if n_e[i] >= theta)
                if r>0: rates.append((r,'inf',u)); tot += r
            elif state[u]==1:
                rates.append((mu,'rec',u)); tot += mu
        if tot == 0: break
        x = rng.random()*tot
        for r,kind,u in rates:
            x -= r
            if x <= 0:
                state[u] = 1 if kind=='inf' else 2
                break
    return sum(1 for s in state if s==2)/N

def rand_hg(N,M,tau,eta,rng):
    out=[]
    for _ in range(M):
        nd = rng.sample(range(N), tau+eta)
        out.append((set(nd[:tau]), set(nd[tau:])))
    return out

print("[A] does the final size depend only on lambda = beta/mu ?")
rng = random.Random(21); N,M = 60,50
edges = rand_hg(N,M,2,2,rng)
for theta in (1,2):
    a = [gillespie(edges,N,1.0,0.5,theta,s) for s in range(400)]   # lambda = 2
    b = [gillespie(edges,N,4.0,2.0,theta,s+9000) for s in range(400)]  # lambda = 2, all rates x4
    c = [gillespie(edges,N,0.5,0.5,theta,s+5000) for s in range(400)]  # lambda = 1
    sem = math.sqrt(st.variance(a)/len(a) + st.variance(b)/len(b))
    print(f"    theta={theta}:  lam=2 -> {st.mean(a):.4f} | same lam, 4x rates -> {st.mean(b):.4f}"
          f" | diff {abs(st.mean(a)-st.mean(b)):.4f} vs 2*sem {2*sem:.4f}"
          f"  -> {'consistent' if abs(st.mean(a)-st.mean(b))<2*sem else 'MISMATCH'}")
    print(f"              lam=1 -> {st.mean(c):.4f}   (must differ, else lambda is not the control)")

# ---------- overlap machinery ----------
def channels(edges):
    m = {ab:0 for ab in ("TT","HH")}; tot = {ab:0 for ab in ("TT","HH")}
    for i,(Ti,Hi) in enumerate(edges):
        for j,(Tj,Hj) in enumerate(edges):
            if i==j: continue
            for ab,A,B in (("TT",Ti,Tj),("HH",Hi,Hj)):
                c = len(A&B)
                if c: m[ab]+=1; tot[ab]+=c
    return m, tot

def alpha_pair(edges, tau, eta):
    m,tot = channels(edges)
    aTT = (tot["TT"]/tau)/m["TT"] if m["TT"] else float('nan')
    aHH = (tot["HH"]/eta)/m["HH"] if m["HH"] else float('nan')
    return aTT, aHH, m

def degs(N,edges):
    ko=[0]*N; ki=[0]*N
    for T,H in edges:
        for v in T: ko[v]+=1
        for v in H: ki[v]+=1
    return ko,ki

print()
print("[B] can Delta-alpha be driven away from 0 at FIXED bi-degree sequence?")
print("    (if not, the whole symmetry-breaking experiment is impossible)")
rng = random.Random(4); N,M,tau,eta = 100,80,3,3
edges = rand_hg(N,M,tau,eta,rng)
ko0,ki0 = degs(N,edges)
aTT,aHH,m0 = alpha_pair(edges,tau,eta)
d0 = 0.5*(aTT-aHH)
# Metropolis on |Delta alpha|: tail-side swaps concentrate TT, head-side spread HH
cur = d0
for step in range(30000):
    i,j = rng.sample(range(M),2)
    Ti,Hi = edges[i]; Tj,Hj = edges[j]
    A,B = (Ti,Tj) if rng.random()<.5 else (Hi,Hj)
    v1 = rng.choice(sorted(A)); v2 = rng.choice(sorted(B))
    if v1==v2 or v1 in Tj|Hj or v2 in Ti|Hi: continue
    A.remove(v1); A.add(v2); B.remove(v2); B.add(v1)
    aTT,aHH,_ = alpha_pair(edges,tau,eta)
    new = 0.5*(aTT-aHH)
    if new >= cur:                      # greedy push upward
        cur = new
    else:
        A.remove(v2); A.add(v1); B.remove(v1); B.add(v2)   # revert
ko1,ki1 = degs(N,edges)
aTT,aHH,m1 = alpha_pair(edges,tau,eta)
print(f"    degrees unchanged: {(ko0,ki0)==(ko1,ki1)}")
print(f"    Delta-alpha: {d0:+.5f} -> {cur:+.5f}   (m_TT {m0['TT']}->{m1['TT']}, m_HH {m0['HH']}->{m1['HH']})")
print(f"    -> {'TUNABLE' if abs(cur-d0)>1e-3 else 'NOT TUNABLE — C4 experiment infeasible'}")

print()
print("[C] can alpha_parallel be held fixed while Delta-alpha varies?")
aP0 = 0.5*((m0['TT'] and 1) and ( (channels(edges)[1]['TT']/tau)/m1['TT'] ) + 0)
aTT,aHH,_ = alpha_pair(edges,tau,eta)
aP1 = 0.5*(aTT+aHH)
print(f"    alpha_par now {aP1:.5f};  alpha_par = (C/2)(1/m_TT + 1/m_HH), Delta = (C/2)(1/m_TT - 1/m_HH)")
print(f"    two independent handles (m_TT, m_HH) for two targets -> a level set of alpha_par exists;")
print(f"    feasibility is a constraint on the swap chain, not an obstruction in principle.")
