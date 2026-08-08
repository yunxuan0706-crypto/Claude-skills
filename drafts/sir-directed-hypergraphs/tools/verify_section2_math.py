#!/usr/bin/env python3
"""Numerical audit of every identity asserted in Section II.

Builds random directed hypergraphs (tails and heads disjoint by
construction) and machine-checks the claims the manuscript makes:

  1  handshake sums                 sum k_out = tau*M,  sum k_in = eta*M
  2  reduction identity             |supp(e) & supp(f)| = N_TT+N_TH+N_HT+N_HH
  3  per-edge sum rules             the four Eq.-(12)-type identities
  4  channel totals                 closed forms in the degree sequence
  5  degree-preserving swaps        totals and alpha*m invariant, alpha moves
  6  reversal parity                r_io even; a^TT<->a^HH; a^HT == a^TH

Exit 0 iff every check passes.
"""
import math
import random
import sys

FAIL = 0


def check(name, ok, detail=""):
    global FAIL
    print(("PASS  " if ok else "FAIL  ") + name + (f"   {detail}" if detail else ""))
    if not ok:
        FAIL += 1


def build(N, M, tau, eta, rng):
    edges = []
    for _ in range(M):
        nodes = rng.sample(range(N), tau + eta)
        edges.append((set(nodes[:tau]), set(nodes[tau:])))
    return edges


def degrees(N, edges):
    ko, ki = [0] * N, [0] * N
    for T, H in edges:
        for v in T:
            ko[v] += 1
        for v in H:
            ki[v] += 1
    return ko, ki


def channel_counts(edges):
    out = {ab: {} for ab in ("TT", "HH", "HT", "TH")}
    for i, (Ti, Hi) in enumerate(edges):
        for j, (Tj, Hj) in enumerate(edges):
            if i == j:
                continue
            for ab, A, B in (("TT", Ti, Tj), ("HH", Hi, Hj), ("HT", Hi, Tj), ("TH", Ti, Hj)):
                c = len(A & B)
                if c:
                    out[ab][(i, j)] = c
    return out


def alpha_stats(nab, ab, tau, eta):
    den = {"TT": tau, "HH": eta, "HT": min(tau, eta), "TH": min(tau, eta)}[ab]
    vals = [c / den for c in nab[ab].values()]
    if not vals:
        return 0.0, 0, 0.0
    return sum(vals) / len(vals), len(vals), sum(vals)


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(x, y))
    vx = sum((a - mx) ** 2 for a in x)
    vy = sum((b - my) ** 2 for b in y)
    return cov / math.sqrt(vx * vy)


def audit(N, M, tau, eta, seed):
    rng = random.Random(seed)
    edges = build(N, M, tau, eta, rng)
    ko, ki = degrees(N, edges)
    tag = f"(N={N},M={M},tau={tau},eta={eta})"

    # 1 -- handshakes
    check(f"handshake {tag}", sum(ko) == tau * M and sum(ki) == eta * M)

    # 2 -- reduction identity on every ordered pair
    ok = True
    for i, (Ti, Hi) in enumerate(edges):
        for j, (Tj, Hj) in enumerate(edges):
            if i == j:
                continue
            lhs = len((Ti | Hi) & (Tj | Hj))
            rhs = (len(Ti & Tj) + len(Ti & Hj) + len(Hi & Tj) + len(Hi & Hj))
            if lhs != rhs:
                ok = False
    check(f"reduction identity {tag}", ok)

    # 3 -- per-edge sum rules
    nab = channel_counts(edges)
    ok = True
    for i, (Ti, Hi) in enumerate(edges):
        sTT = sum(c for (a, b), c in nab["TT"].items() if a == i)
        sHH = sum(c for (a, b), c in nab["HH"].items() if a == i)
        sHT = sum(c for (a, b), c in nab["HT"].items() if a == i)
        sTH = sum(c for (a, b), c in nab["TH"].items() if a == i)
        if sTT != sum(ko[v] - 1 for v in Ti):
            ok = False
        if sHH != sum(ki[v] - 1 for v in Hi):
            ok = False
        if sHT != sum(ko[v] for v in Hi):   # uses T(e) disjoint from H(e)
            ok = False
        if sTH != sum(ki[v] for v in Ti):
            ok = False
    check(f"per-edge sum rules {tag}", ok)

    # 4 -- channel totals against closed forms in the degree sequence
    tot = {ab: sum(nab[ab].values()) for ab in nab}
    check(f"total TT {tag}", tot["TT"] == sum(k * (k - 1) for k in ko))
    check(f"total HH {tag}", tot["HH"] == sum(k * (k - 1) for k in ki))
    check(f"totals HT=TH=sum(kin*kout) {tag}",
          tot["HT"] == tot["TH"] == sum(a * b for a, b in zip(ki, ko)))

    # 5 -- degree-preserving double swaps: totals and alpha*m invariant
    a0 = {ab: alpha_stats(nab, ab, tau, eta) for ab in nab}
    swapped = 0
    for _ in range(4000):
        if swapped >= 300:
            break
        i, j = rng.sample(range(M), 2)
        Ti, Hi = edges[i]
        Tj, Hj = edges[j]
        side = rng.random() < 0.5
        A, B = (Ti, Tj) if side else (Hi, Hj)
        v1 = rng.choice(sorted(A))
        v2 = rng.choice(sorted(B))
        if v1 == v2 or v1 in Tj | Hj or v2 in Ti | Hi:
            continue
        A.remove(v1); A.add(v2)
        B.remove(v2); B.add(v1)
        swapped += 1
    ko2, ki2 = degrees(N, edges)
    check(f"swaps preserve degrees {tag}", (ko2, ki2) == (ko, ki), f"{swapped} swaps")
    nab2 = channel_counts(edges)
    a1 = {ab: alpha_stats(nab2, ab, tau, eta) for ab in nab2}
    ok_tot = all(abs(a0[ab][2] - a1[ab][2]) < 1e-9 for ab in a0)
    check(f"sum of n_ab (=alpha*m) invariant under swaps {tag}", ok_tot,
          "C_ab fixed by degree sequence")
    moved = max(abs(a0[ab][0] - a1[ab][0]) for ab in a0)
    check(f"alpha itself moved via m {tag}", moved > 1e-4, f"max |d alpha|={moved:.4f}")

    # 6 -- reversal parity
    edgesR = [(H, T) for (T, H) in edges]
    koR, kiR = degrees(N, edgesR)
    check(f"reversal swaps degree roles {tag}", koR == ki2 and kiR == ko2)
    check(f"r_io reversal-even {tag}",
          abs(pearson(ki2, ko2) - pearson(kiR, koR)) < 1e-12)
    nabR = channel_counts(edgesR)
    aR = {ab: alpha_stats(nabR, ab, eta, tau) for ab in nabR}
    check(f"alpha^TT(R H)=alpha^HH(H) {tag}", abs(aR["TT"][0] - a1["HH"][0]) < 1e-12)
    check(f"alpha^HH(R H)=alpha^TT(H) {tag}", abs(aR["HH"][0] - a1["TT"][0]) < 1e-12)
    check(f"alpha^HT == alpha^TH on the same graph {tag}",
          abs(a1["HT"][0] - a1["TH"][0]) < 1e-12,
          "serial channels coincide pair-by-relabelling")


if __name__ == "__main__":
    audit(N=120, M=90, tau=3, eta=2, seed=7)
    audit(N=150, M=120, tau=2, eta=2, seed=11)
    print()
    print("ALL CHECKS PASSED" if FAIL == 0 else f"{FAIL} CHECK(S) FAILED")
    sys.exit(1 if FAIL else 0)
