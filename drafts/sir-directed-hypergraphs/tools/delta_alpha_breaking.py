#!/usr/bin/env python3
"""Decisive test: does Delta-alpha alone break direction symmetry?

The earlier pilot was underpowered by design -- one structure, N=160, deep
in the supercritical regime where Pi saturates and its sensitivity to
structure is worst. This version fixes all three:

  * runs NEAR THRESHOLD, where Pi is most sensitive;
  * larger N;
  * averages over K independent structures, so the reported error bar is an
    ensemble error bar rather than one structure's fluctuation.

Design. Build a bi-degree sequence that is swap-symmetric by construction
(mirror pairing), so the sequence contributes no reversal-odd part and
Delta-alpha is the only candidate source. Drive Delta-alpha up with
degree-preserving swaps, then compare Pi and S_bar on H against R(H) -- the
same structure with Delta-alpha of the opposite sign. alpha_parallel is
reversal-even, so its drift cannot fake a paired difference.

Overlap counts are maintained incrementally: a double swap touches two
hyperedges, so only their rows need recomputing, making the drive loop
O(M) per proposal instead of O(M^2).
"""
import math
import random
import statistics as st
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from feasibility_check import index, gillespie_fast, joint_of, mf_threshold


def mirror_structure(Nh, Mh, tau, eta, rng):
    """Hypergraph on 2*Nh nodes whose bi-degree sequence is swap-symmetric."""
    half = []
    for _ in range(Mh):
        nd = rng.sample(range(Nh), tau + eta)
        half.append((set(nd[:tau]), set(nd[tau:])))
    return half + [({v + Nh for v in H}, {v + Nh for v in T}) for T, H in half]


class Overlap:
    """Incrementally maintained TT / HH channel counts."""

    def __init__(self, edges, tau, eta):
        self.e, self.tau, self.eta = edges, tau, eta
        self.TT, self.HH = {}, {}
        for i in range(len(edges)):
            self.row(i)

    def row(self, i):
        Ti, Hi = self.e[i]
        for j in range(len(self.e)):
            if j == i:
                continue
            Tj, Hj = self.e[j]
            c = len(Ti & Tj)
            if c:
                self.TT[(i, j)] = self.TT[(j, i)] = c
            else:
                self.TT.pop((i, j), None)
                self.TT.pop((j, i), None)
            c = len(Hi & Hj)
            if c:
                self.HH[(i, j)] = self.HH[(j, i)] = c
            else:
                self.HH.pop((i, j), None)
                self.HH.pop((j, i), None)

    def alphas(self):
        aTT = (sum(self.TT.values()) / self.tau) / len(self.TT)
        aHH = (sum(self.HH.values()) / self.eta) / len(self.HH)
        return aTT, aHH

    def delta(self):
        a, b = self.alphas()
        return 0.5 * (a - b)

    def par(self):
        a, b = self.alphas()
        return 0.5 * (a + b)


def drive(edges, ov, rng, proposals, target=None):
    """Greedy degree-preserving swaps that increase Delta-alpha."""
    M = len(edges)
    best = ov.delta()
    for _ in range(proposals):
        if target is not None and best >= target:
            break
        i, j = rng.sample(range(M), 2)
        Ti, Hi = edges[i]
        Tj, Hj = edges[j]
        A, B = (Ti, Tj) if rng.random() < 0.5 else (Hi, Hj)
        v1 = rng.choice(sorted(A))
        v2 = rng.choice(sorted(B))
        if v1 == v2 or v1 in Tj | Hj or v2 in Ti | Hi:
            continue
        A.remove(v1); A.add(v2)
        B.remove(v2); B.add(v1)
        ov.row(i); ov.row(j)
        d = ov.delta()
        if d > best:
            best = d
        else:
            A.remove(v2); A.add(v1)
            B.remove(v1); B.add(v2)
            ov.row(i); ov.row(j)
    return best


def sizes(N, edges, lam, runs, seed0):
    tails, heads, tails_of = index(N, edges)
    return [gillespie_fast(N, tails, heads, tails_of, lam, 1,
                           random.Random(seed0 + s),
                           [random.Random(seed0 + 7919 * s).randrange(N)])
            for s in range(runs)]


def summarise(S, cut):
    big = [x for x in S if x > cut]
    return len(big) / len(S), (st.mean(big) if big else float("nan"))


def main(N_half=200, M_half=150, tau=2, eta=2, K=6, runs=5000,
         ratios=(1.0, 1.3, 1.7), target=0.07, seed=404):
    rng = random.Random(seed)
    N = 2 * N_half
    print(f"N={N}  M={2*M_half}  tau=eta={tau}  structures={K}  runs/point={runs}")

    # locate the threshold once, on a representative structure
    probe = mirror_structure(N_half, M_half, tau, eta, random.Random(seed))
    jc, _, _ = joint_of(N, probe)
    lc = mf_threshold(jc, tau)
    print(f"lambda_c(mean-field reference) = {lc:.4f}\n")

    pooled = {r: [] for r in ratios}
    for k in range(K):
        edges = mirror_structure(N_half, M_half, tau, eta, random.Random(seed + 31 * k))
        ov = Overlap(edges, tau, eta)
        d0 = ov.delta()
        d1 = drive(edges, ov, random.Random(seed + 97 * k), 40000, target)
        edgesR = [(H, T) for T, H in edges]
        print(f"structure {k}: Delta-alpha {d0:+.4f} -> {d1:+.4f}   par={ov.par():.4f}")
        for r in ratios:
            lam = r * lc
            SA = sizes(N, edges, lam, runs, 100000 + 1000 * k)
            SB = sizes(N, edgesR, lam, runs, 100000 + 1000 * k)
            cut = 0.02
            pA, sA = summarise(SA, cut)
            pB, sB = summarise(SB, cut)
            sep = math.sqrt(pA * (1 - pA) / runs + pB * (1 - pB) / runs)
            pooled[r].append((pA - pB, sep, pA, pB, sA, sB,
                              st.mean(SA), st.mean(SB)))
            print(f"   lam={r:.2f}*lc: Pi {pA:.4f} vs {pB:.4f}"
                  f"  dPi={pA-pB:+.4f}+/-{sep:.4f}  ({abs(pA-pB)/sep:.2f} sig)"
                  f"   <S> {st.mean(SA):.4f}/{st.mean(SB):.4f}")

    print("\n" + "=" * 66)
    print("POOLED OVER STRUCTURES (ensemble error bar)")
    print("=" * 66)
    for r in ratios:
        d = [x[0] for x in pooled[r]]
        md = st.mean(d)
        sem = math.sqrt(st.variance(d) / len(d)) if len(d) > 1 else float("nan")
        mS = st.mean([x[6] - x[7] for x in pooled[r]])
        print(f"lam={r:.2f}*lc:  <dPi> = {md:+.5f} +/- {sem:.5f}"
              f"   ({abs(md)/sem:.2f} sigma)   <d<S>> = {mS:+.5f}"
              f"   -> {'BREAKS' if abs(md) > 3*sem else 'not resolved'}")


if __name__ == "__main__":
    main()
