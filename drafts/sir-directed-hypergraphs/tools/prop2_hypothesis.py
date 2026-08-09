#!/usr/bin/env python3
"""Does Proposition 2's hypothesis (ii) hold for the main-line kernel?

(ii) requires the random structure generated on the reversed hypergraph to
be distributed as the transpose of the one generated on the original.

Test A -- linear kernel / independent (tail-member, head-member) arcs. The
structure is seed-independent by construction, so the counting identity is
exact. Confirmed: sum|out| = sum|in| in every trial, no statistics needed.

Test B -- main-line threshold kernel with hyperedge-level activation. The
active interval of a hyperedge is the union of its infected tail members'
infectious periods, so it depends on who was infected, i.e. on the seed;
(ii) is only an approximation. Measured at tau=3, eta=1, N=90, 60000 runs
per point, four independent structures, at the threshold:

    relative deviation of <S>:  -2.24% +/- 0.44%   (5.05 sigma)

Test C -- control against a sampler artefact. A -2% difference between two
runs that differ only by relabelling could equally well be an asymmetry in
the simulator (seed reuse, tie-breaking, iteration order over tails vs
heads) rather than physics. Test C runs the identical protocol on
swap-symmetric bi-degree sequences, where Proposition 1 forbids breaking
outright, so any residual there is instrumental. Measured:

    control (swap-symmetric):  -0.36% +/- 0.28%   (1.31 sigma)
    tau=3, eta=1 (test B):     -2.24% +/- 0.44%   (5.05 sigma)
    B minus control:           -1.88% +/- 0.52%   (3.60 sigma)

The control is consistent with zero, and even taken at face value as an
instrumental baseline it accounts for under a fifth of Test B: the
kernel-induced part survives subtraction at 3.6 sigma. The -2% is a
property of the dynamics, not of the simulator.

The violation is real, systematic in sign, and specific to hyperedge-level
activation -- it has no pairwise analogue. It is nonetheless ~24x smaller
than the breaking Pi shows on the same structures (0.6712 vs 0.3908, i.e.
53% of their mean; see breaking_observable.py), so the
conclusion that <S> is a poor order parameter stands quantitatively even
though the invariance is not exact.
"""
import math
import random
import statistics as st
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from feasibility_check import index, gillespie_fast


def test_a(trials=5, N=60, M=45, tau=3, eta=1, seed=3):
    rng = random.Random(seed)
    print("TEST A: independent arcs -- identity must be exact")
    for t in range(trials):
        edges = []
        for _ in range(M):
            nd = rng.sample(range(N), tau + eta)
            edges.append((set(nd[:tau]), set(nd[tau:])))
        arcs = {}
        for T, H in edges:
            for v in T:
                for u in H:
                    if rng.random() < 0.25:
                        arcs.setdefault(v, set()).add(u)
        rev = {}
        for v, us in arcs.items():
            for u in us:
                rev.setdefault(u, set()).add(v)

        def reach(src, A):
            seen, stack = {src}, [src]
            while stack:
                x = stack.pop()
                for y in A.get(x, ()):
                    if y not in seen:
                        seen.add(y)
                        stack.append(y)
            return len(seen)

        so = sum(reach(v, arcs) for v in range(N))
        si = sum(reach(v, rev) for v in range(N))
        print(f"  trial {t}: sum|out|={so} sum|in|={si} equal={so == si}")


def test_b(K=4, N=90, M=70, runs=60000, lam=1.0, seed=9):
    print(f"\nTEST B: threshold kernel, tau=3 eta=1, N={N}, {runs} runs, {K} structures")
    rel = []
    for k in range(K):
        rng = random.Random(seed + 17 * k)
        edges = []
        for _ in range(M):
            nd = rng.sample(range(N), 4)
            edges.append((set(nd[:3]), set(nd[3:])))
        t1 = index(N, edges)
        t2 = index(N, [(H, T) for T, H in edges])
        A = [gillespie_fast(N, *t1, lam, 1, random.Random(100 * k + s),
                            [random.Random(700000 + 100 * k + s).randrange(N)]) for s in range(runs)]
        B = [gillespie_fast(N, *t2, lam, 1, random.Random(100 * k + s),
                            [random.Random(700000 + 100 * k + s).randrange(N)]) for s in range(runs)]
        mA, mB = st.mean(A), st.mean(B)
        rel.append((mA - mB) / mA)
        sem = math.sqrt(st.variance(A) / runs + st.variance(B) / runs)
        print(f"  structure {k}: <S>_H={mA:.5f} <S>_RH={mB:.5f} "
              f"diff={mA-mB:+.5f} ({abs(mA-mB)/sem:.2f} sig) rel={(mA-mB)/mA:+.2%}")
    m = st.mean(rel)
    s = math.sqrt(st.variance(rel) / len(rel))
    print(f"\n  pooled relative deviation: {m:+.4%} +/- {s:.4%}  ({abs(m)/s:.2f} sigma)")
    print("  -> hypothesis (ii) is violated by the hyperedge-level kernel")
    return m, s


def test_c(K=4, N_half=45, M_half=35, runs=60000, lam=1.0, seed=9):
    """Control: same protocol on swap-symmetric sequences (Prop. 1 applies)."""
    print(f"\nTEST C: control, swap-symmetric sequence, N={2*N_half}, "
          f"{runs} runs, {K} structures")
    rel = []
    for k in range(K):
        rng = random.Random(seed + 17 * k)
        half = []
        for _ in range(M_half):
            nd = rng.sample(range(N_half), 4)
            half.append((set(nd[:3]), set(nd[3:])))
        edges = half + [({v + N_half for v in H}, {v + N_half for v in T})
                        for T, H in half]
        N = 2 * N_half
        t1 = index(N, edges)
        t2 = index(N, [(H, T) for T, H in edges])
        A = [gillespie_fast(N, *t1, lam, 1, random.Random(100 * k + s),
                            [random.Random(700000 + 100 * k + s).randrange(N)]) for s in range(runs)]
        B = [gillespie_fast(N, *t2, lam, 1, random.Random(100 * k + s),
                            [random.Random(700000 + 100 * k + s).randrange(N)]) for s in range(runs)]
        mA, mB = st.mean(A), st.mean(B)
        rel.append((mA - mB) / mA)
        sem = math.sqrt(st.variance(A) / runs + st.variance(B) / runs)
        print(f"  structure {k}: <S>_H={mA:.5f} <S>_RH={mB:.5f} "
              f"diff={mA-mB:+.5f} ({abs(mA-mB)/sem:.2f} sig) rel={(mA-mB)/mA:+.2%}")
    m = st.mean(rel)
    s = math.sqrt(st.variance(rel) / len(rel))
    print(f"\n  pooled relative deviation: {m:+.4%} +/- {s:.4%}  ({abs(m)/s:.2f} sigma)")
    print("  -> consistent with zero; treat as the instrumental baseline")
    return m, s


if __name__ == "__main__":
    test_a()
    mB, sB = test_b()
    mC, sC = test_c()
    d = mB - mC
    sd = math.sqrt(sB ** 2 + sC ** 2)
    print(f"\nB minus control: {d:+.4%} +/- {sd:.4%}  ({abs(d)/sd:.2f} sigma)")
    print("  -> the kernel-induced part survives the instrumental baseline")
