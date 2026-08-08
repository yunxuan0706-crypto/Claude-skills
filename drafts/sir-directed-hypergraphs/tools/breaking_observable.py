"""Which observable actually detects direction-symmetry breaking?

Result (see Section II.D, Propositions 1 and 2):

  * the MEAN outbreak size from a uniform seed does NOT break, on any
    structure -- sum_v |out(v)| = sum_v |in(v)| counts the same ordered
    reachable pairs, and reversal only transposes them. Measured: <= 1.4
    sigma even at tau=3, eta=1, and even with Delta-alpha driven to 0.114;
  * the epidemic PROBABILITY and the CONDITIONAL final size break hard and
    in opposite directions, their product staying pinned to the invariant
    mean. Measured at tau=3, eta=1, lambda=2.4:
        <S>       0.2928 vs 0.2883   (0.57 sigma -- no break)
        P(large)  0.6712 vs 0.3908   (20.7 sigma -- breaks)
        <S|large> 0.4230 vs 0.7065   (breaks the other way)
        product   0.2839 vs 0.2761   (compensation holds)

So Section VI must use (Pi, S_bar), not <S>, and should report <S> as an
internal consistency check.
"""
import math, random, statistics as st

def rand_hg_symmetric(N, M, tau, eta, rng):
    """Random hypergraph whose bi-degree sequence is swap-symmetric by
    construction: build it, then pair each node with a partner carrying the
    mirrored (k_in,k_out). Simplest route: build any hypergraph, then use its
    reversal's degree sequence too -- i.e. take the disjoint union of H and R(H)."""
    edges = []
    for _ in range(M):
        nd = rng.sample(range(N // 2), tau + eta)
        edges.append((set(nd[:tau]), set(nd[tau:])))
    # mirror copy on the other half of the node set, with tails/heads swapped
    off = N // 2
    mirror = [({v + off for v in H}, {v + off for v in T}) for T, H in edges]
    return edges + mirror

def degrees(N, edges):
    ko, ki = [0]*N, [0]*N
    for T, H in edges:
        for v in T: ko[v] += 1
        for v in H: ki[v] += 1
    return ko, ki

def alphas(edges, tau, eta):
    mTT = mHH = 0; tTT = tHH = 0
    for i, (Ti, Hi) in enumerate(edges):
        for j, (Tj, Hj) in enumerate(edges):
            if i == j: continue
            c = len(Ti & Tj)
            if c: mTT += 1; tTT += c
            c = len(Hi & Hj)
            if c: mHH += 1; tHH += c
    aTT = (tTT/tau)/mTT if mTT else 0.0
    aHH = (tHH/eta)/mHH if mHH else 0.0
    return aTT, aHH, mTT, mHH

def gillespie(edges, N, heads, lam, theta, seed, nseed):
    rng = random.Random(seed)
    state = bytearray(N)
    inf = set(rng.sample(range(N), nseed))
    for v in inf: state[v] = 1
    n_e = [sum(1 for v in T if state[v] == 1) for T, _ in edges]
    while inf:
        # rates: infection lam per active incident hyperedge, recovery 1
        cand = []
        tot = 0.0
        touched = set()
        for i, (T, H) in enumerate(edges):
            if n_e[i] >= theta:
                for u in H:
                    if state[u] == 0: touched.add(u)
        for u in touched:
            r = lam*sum(1 for i in heads[u] if n_e[i] >= theta)
            cand.append((r, 1, u)); tot += r
        for u in inf:
            cand.append((1.0, 2, u)); tot += 1.0
        x = rng.random()*tot
        for r, kind, u in cand:
            x -= r
            if x <= 0:
                if kind == 1:
                    state[u] = 1; inf.add(u)
                    for i, (T, _) in enumerate(edges):
                        if u in T: n_e[i] += 1
                else:
                    state[u] = 2; inf.discard(u)
                    for i, (T, _) in enumerate(edges):
                        if u in T: n_e[i] -= 1
                break
    return sum(1 for s in state if s == 2)/N

def head_index(N, edges):
    h = [[] for _ in range(N)]
    for i, (_, H) in enumerate(edges):
        for v in H: h[v].append(i)
    return h

def mean_S(edges, N, lam, theta, nseed, runs, seed0):
    h = head_index(N, edges)
    return [gillespie(edges, N, h, lam, theta, seed0+s, nseed) for s in range(runs)]

# ---------------------------------------------------------------- run ----
rng = random.Random(101)
N, M, tau, eta = 80, 60, 2, 2
edges = rand_hg_symmetric(N, M, tau, eta, rng)
M = len(edges)
ko, ki = degrees(N, edges)
seq = sorted(zip(ki, ko))
print(f"N={N}  M={M}  tau=eta={tau}")
print(f"bi-degree sequence swap-symmetric: {seq == sorted((b,a) for a,b in seq)}")
aTT, aHH, mTT, mHH = alphas(edges, tau, eta)
print(f"start:  Delta-alpha = {0.5*(aTT-aHH):+.5f}   alpha_par = {0.5*(aTT+aHH):.5f}")

# drive Delta-alpha up while keeping alpha_parallel near its start value
target_par = 0.5*(aTT+aHH)
best = 0.5*(aTT-aHH)
for step in range(60000):
    i, j = rng.sample(range(M), 2)
    Ti, Hi = edges[i]; Tj, Hj = edges[j]
    A, B = (Ti, Tj) if rng.random() < .5 else (Hi, Hj)
    v1 = rng.choice(sorted(A)); v2 = rng.choice(sorted(B))
    if v1 == v2 or v1 in Tj | Hj or v2 in Ti | Hi: continue
    A.remove(v1); A.add(v2); B.remove(v2); B.add(v1)
    aTT, aHH, _, _ = alphas(edges, tau, eta)
    d = 0.5*(aTT-aHH); p = 0.5*(aTT+aHH)
    if d > best and abs(p - target_par) < 0.02:
        best = d
    else:
        A.remove(v2); A.add(v1); B.remove(v1); B.add(v2)
ko2, ki2 = degrees(N, edges)
aTT, aHH, mTT, mHH = alphas(edges, tau, eta)
print(f"after:  Delta-alpha = {0.5*(aTT-aHH):+.5f}   alpha_par = {0.5*(aTT+aHH):.5f}"
      f"   (degrees unchanged: {(ko2,ki2)==(ko,ki)})")

edgesR = [(H, T) for T, H in edges]
aTTr, aHHr, _, _ = alphas(edgesR, tau, eta)
print(f"reversed: Delta-alpha = {0.5*(aTTr-aHHr):+.5f}  (must be the negative)")

print()
for theta, nseed, lam in ((1, 1, 1.2),):
    RUNS = 1500
    a = mean_S(edges,  N, lam, theta, nseed, RUNS, 1000)
    b = mean_S(edgesR, N, lam, theta, nseed, RUNS, 1000)   # same seeds
    diff = [x-y for x, y in zip(a, b)]
    md = st.mean(diff); sem = math.sqrt(st.variance(diff)/len(diff))
    print(f"theta={theta}, lambda={lam}, seeds={nseed}, runs={RUNS}")
    print(f"  <S> on H    = {st.mean(a):.5f}")
    print(f"  <S> on R(H) = {st.mean(b):.5f}")
    print(f"  paired diff = {md:+.5f} +/- {sem:.5f}  ({abs(md)/sem:.2f} sigma)")
    print(f"  -> {'BREAKING DETECTED' if abs(md) > 3*sem else 'no breaking at this size/parameters'}")
