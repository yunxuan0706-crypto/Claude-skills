#!/usr/bin/env python3
"""Pre-run of the execution plan's remaining feasibility questions.

  W1  Is there room for C1? Compare the hyperdegree mean-field final size
      (annealed, Li-et-al style classes) against exact Gillespie near the
      threshold, at two system sizes. If mean-field were already accurate,
      the plan's early-stop would fire before EBCM is even derived.
  W4  Is S_cut well defined? Look for a gap in the outbreak-size histogram.
  W5  mu->0 endpoint: at large lambda the final infected set should be the
      deterministic reachability set (theta=1), which is what both limits
      of the C5 bridge describe.
  W4' How far can Delta-alpha be pushed while alpha_parallel is pinned?

Also contains the incremental-rate Gillespie sampler week 1 needs; it is
validated statistically against the simple O(M)-per-event reference
implementation before being used.
"""
import math
import random
import statistics as st

# ------------------------------------------------------------ builders ----
def build(N, M, tau, eta, rng):
    edges = []
    for _ in range(M):
        nd = rng.sample(range(N), tau + eta)
        edges.append((set(nd[:tau]), set(nd[tau:])))
    return edges


def index(N, edges):
    tails_of = [[] for _ in range(N)]
    heads = []
    tails = []
    for i, (T, H) in enumerate(edges):
        tails.append(list(T))
        heads.append(list(H))
        for v in T:
            tails_of[v].append(i)
    return tails, heads, tails_of


# --------------------------------------------------- fast exact Gillespie ----
def gillespie_fast(N, tails, heads, tails_of, lam, theta, rng, seeds):
    state = bytearray(N)
    n_e = [0] * len(tails)
    pressure = {}
    sump = 0
    inf_list, inf_pos = [], {}

    def activate(ei):
        nonlocal sump
        for w in heads[ei]:
            if state[w] == 0:
                pressure[w] = pressure.get(w, 0) + 1
                sump += 1

    def deactivate(ei):
        nonlocal sump
        for w in heads[ei]:
            if state[w] == 0:
                c = pressure[w] - 1
                if c:
                    pressure[w] = c
                else:
                    del pressure[w]
                sump -= 1

    def infect(u):
        nonlocal sump
        if u in pressure:
            sump -= pressure.pop(u)
        state[u] = 1
        inf_pos[u] = len(inf_list)
        inf_list.append(u)
        for ei in tails_of[u]:
            n_e[ei] += 1
            if n_e[ei] == theta:
                activate(ei)

    def recover(u):
        state[u] = 2
        i = inf_pos.pop(u)
        last = inf_list.pop()
        if last != u:
            inf_list[i] = last
            inf_pos[last] = i
        for ei in tails_of[u]:
            if n_e[ei] == theta:
                deactivate(ei)
            n_e[ei] -= 1

    for v in seeds:
        infect(v)
    while inf_list:
        tot = lam * sump + len(inf_list)
        if rng.random() * tot < lam * sump:
            x = rng.random() * sump
            for u, c in pressure.items():
                x -= c
                if x <= 0:
                    infect(u)
                    break
        else:
            recover(rng.choice(inf_list))
    return sum(1 for s in state if s == 2) / N


def gillespie_ref(edges, N, lam, theta, rng, seeds):
    """Simple reference implementation (recomputes everything each event)."""
    state = [0] * N
    for v in seeds:
        state[v] = 1
    heads_of = [[] for _ in range(N)]
    for i, (T, H) in enumerate(edges):
        for v in H:
            heads_of[v].append(i)
    while True:
        n_e = [sum(1 for v in T if state[v] == 1) for T, _ in edges]
        cand, tot = [], 0.0
        for u in range(N):
            if state[u] == 0:
                r = lam * sum(1 for i in heads_of[u] if n_e[i] >= theta)
                if r > 0:
                    cand.append((r, 1, u))
                    tot += r
            elif state[u] == 1:
                cand.append((1.0, 2, u))
                tot += 1.0
        if tot == 0:
            break
        x = rng.random() * tot
        for r, kind, u in cand:
            x -= r
            if x <= 0:
                state[u] = kind
                break
    return sum(1 for s in state if s == 2) / N


# --------------------------------------------- hyperdegree mean field ----
def meanfield_final(joint_classes, tau, lam, eps=1e-4, dt=0.02, tmax=400.0):
    """Annealed hyperdegree mean field for theta=1 (Li-et-al style classes)."""
    cls = [(ki, ko, P) for (ki, ko), P in joint_classes.items()]
    kom = sum(ko * P for _, ko, P in cls)
    s = [1.0 - eps] * len(cls)
    i = [eps] * len(cls)
    r = [0.0] * len(cls)
    t = 0.0
    while t < tmax:
        pT = sum(ko * P * iv for (ki, ko, P), iv in zip(cls, i)) / kom
        a = 1.0 - (1.0 - pT) ** tau
        itot = 0.0
        for idx, (ki, ko, P) in enumerate(cls):
            inf_rate = lam * ki * a * s[idx]
            s[idx] -= inf_rate * dt
            r[idx] += i[idx] * dt
            i[idx] += (inf_rate - i[idx]) * dt
            itot += i[idx] * P
        t += dt
        if itot < 1e-10:
            break
    return sum(P * (rv + iv) for (ki, ko, P), rv, iv in zip(cls, r, i))


def mf_threshold(joint_classes, tau):
    kom = sum(ko * P for (ki, ko), P in joint_classes.items())
    kk = sum(ki * ko * P for (ki, ko), P in joint_classes.items())
    return kom / (tau * kk)


# ------------------------------------------------------- experiments ----
def joint_of(N, edges):
    from collections import Counter
    ko, ki = [0] * N, [0] * N
    for T, H in edges:
        for v in T:
            ko[v] += 1
        for v in H:
            ki[v] += 1
    c = Counter(zip(ki, ko))
    return {k: v / N for k, v in c.items()}, ko, ki


def w1_meanfield_room():
    print("=" * 72)
    print("W1  mean-field error near threshold (the plan's top risk, pre-run)")
    print("=" * 72)
    tau = eta = 2
    theta = 1
    CUT = 0.05
    for N, runs in ((150, 400), (300, 250)):
        M = int(0.75 * N)
        rng = random.Random(37 + N)
        edges = build(N, M, tau, eta, rng)
        tails, heads, tails_of = index(N, edges)
        jc, _, _ = joint_of(N, edges)
        lc = mf_threshold(jc, tau)
        print(f"\nN={N} M={M}: lambda_c(MF) = {lc:.4f}")
        print(f"{'lam/lc':>7} {'Pi':>7} {'S|large':>9} {'S_MF':>8} {'rel.err':>9}")
        for f in (0.8, 1.0, 1.25, 1.6, 2.2):
            lam = f * lc
            S = [gillespie_fast(N, tails, heads, tails_of, lam, theta,
                                random.Random(1000 + s), [random.Random(500 + s).randrange(N)])
                 for s in range(runs)]
            big = [x for x in S if x > CUT]
            Pi = len(big) / len(S)
            Sb = st.mean(big) if big else float("nan")
            Smf = meanfield_final(jc, tau, lam)
            rel = (Smf - Sb) / Sb if big else float("inf")
            print(f"{f:7.2f} {Pi:7.3f} {Sb:9.4f} {Smf:8.4f} {rel:+9.1%}")


def w4_bimodality():
    print()
    print("=" * 72)
    print("W4  is S_cut well defined? (mass in the would-be gap 0.05-0.15)")
    print("=" * 72)
    tau = eta = 2
    N, M = 300, 225
    rng = random.Random(11)
    edges = build(N, M, tau, eta, rng)
    tails, heads, tails_of = index(N, edges)
    jc, _, _ = joint_of(N, edges)
    lc = mf_threshold(jc, tau)
    for f in (1.25, 1.6, 2.2):
        lam = f * lc
        S = [gillespie_fast(N, tails, heads, tails_of, lam, 1,
                            random.Random(3000 + s), [random.Random(700 + s).randrange(N)])
             for s in range(400)]
        gap = sum(1 for x in S if 0.05 < x < 0.15) / len(S)
        print(f"  lam={f:.2f}*lc: gap mass {gap:.3f}  small {sum(1 for x in S if x<=0.05)/len(S):.3f}"
              f"  large {sum(1 for x in S if x>=0.15)/len(S):.3f}")


def w5_mu0_endpoint():
    print()
    print("=" * 72)
    print("W5  mu->0 endpoint: Gillespie at lam=500 vs deterministic reachability")
    print("=" * 72)
    N, M, tau, eta = 120, 90, 2, 2
    rng = random.Random(13)
    edges = build(N, M, tau, eta, rng)
    tails, heads, tails_of = index(N, edges)

    def reach(v):
        infected = {v}
        fired = [False] * M
        stack = [v]
        while stack:
            u = stack.pop()
            for ei in tails_of[u]:
                if not fired[ei]:
                    fired[ei] = True
                    for w in heads[ei]:
                        if w not in infected:
                            infected.add(w)
                            stack.append(w)
        return len(infected) / N

    diffs = []
    for s in range(150):
        v = random.Random(40 + s).randrange(N)
        g = gillespie_fast(N, tails, heads, tails_of, 500.0, 1, random.Random(9000 + s), [v])
        diffs.append(g - reach(v))
    print(f"  mean |S_gillespie - S_reach| = {st.mean(map(abs,diffs)):.4f}"
          f"   max = {max(map(abs,diffs)):.4f}   (prediction: ~1/(1+lam) per link)")


def w4_constrained_delta():
    print()
    print("=" * 72)
    print("W4' Delta-alpha range with alpha_parallel pinned (open item)")
    print("=" * 72)
    N, M0, tau, eta = 100, 60, 3, 3
    rng = random.Random(101)
    half = []
    for _ in range(M0):
        nd = rng.sample(range(N // 2), tau + eta)
        half.append((set(nd[:tau]), set(nd[tau:])))
    off = N // 2
    edges = half + [({v + off for v in H}, {v + off for v in T}) for T, H in half]
    M = len(edges)

    def alphas():
        mTT = mHH = tTT = tHH = 0
        for i, (Ti, Hi) in enumerate(edges):
            for j, (Tj, Hj) in enumerate(edges):
                if i == j:
                    continue
                c = len(Ti & Tj)
                if c:
                    mTT += 1
                    tTT += c
                c = len(Hi & Hj)
                if c:
                    mHH += 1
                    tHH += c
        return (tTT / tau) / mTT, (tHH / eta) / mHH

    aTT, aHH = alphas()
    par0, d = 0.5 * (aTT + aHH), 0.5 * (aTT - aHH)
    print(f"  start: Delta={d:+.4f}  par={par0:.4f}")
    best = d
    for _ in range(40000):
        i, j = rng.sample(range(M), 2)
        Ti, Hi = edges[i]
        Tj, Hj = edges[j]
        A, B = (Ti, Tj) if rng.random() < .5 else (Hi, Hj)
        v1 = rng.choice(sorted(A))
        v2 = rng.choice(sorted(B))
        if v1 == v2 or v1 in Tj | Hj or v2 in Ti | Hi:
            continue
        A.remove(v1); A.add(v2)
        B.remove(v2); B.add(v1)
        aTT, aHH = alphas()
        nd_, np_ = 0.5 * (aTT - aHH), 0.5 * (aTT + aHH)
        if nd_ > best and abs(np_ - par0) <= 0.02:
            best = nd_
        else:
            A.remove(v2); A.add(v1)
            B.remove(v1); B.add(v2)
    print(f"  best Delta-alpha with |par-par0|<=0.02:  {best:+.4f}"
          f"   (unconstrained benchmark at tau=eta=2 was +0.114)")


if __name__ == "__main__":
    # validate the fast sampler against the reference first
    N, M, tau, eta = 80, 60, 2, 2
    rng = random.Random(5)
    edges = build(N, M, tau, eta, rng)
    tails, heads, tails_of = index(N, edges)
    a = [gillespie_fast(N, tails, heads, tails_of, 1.0, 1,
                        random.Random(s), [random.Random(70 + s).randrange(N)]) for s in range(400)]
    b = [gillespie_ref(edges, N, 1.0, 1,
                       random.Random(50000 + s), [random.Random(70 + s).randrange(N)]) for s in range(400)]
    sem = math.sqrt(st.variance(a) / len(a) + st.variance(b) / len(b))
    print(f"sampler validation: fast {st.mean(a):.4f} vs reference {st.mean(b):.4f}"
          f"  ({abs(st.mean(a)-st.mean(b))/sem:.2f} sigma)  "
          f"{'OK' if abs(st.mean(a)-st.mean(b)) < 3*sem else 'MISMATCH'}")

    w1_meanfield_room()
    w4_bimodality()
    w5_mu0_endpoint()
    w4_constrained_delta()
