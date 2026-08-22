"""Inter-layer group overlap o_ab: construction, measurement, and its effect
on the threshold.

The next-generation matrix (2.7) reads ONLY the joint layer-degree distribution
P(k) and the group sizes m. It is blind, by construction, to how the groups of
different layers sit relative to one another -- to whether a layer-2 group
re-uses a pair of nodes that already share a layer-1 group. That blindness is a
consequence of the local-tree assumption: a shared pair closes a 4-cycle, which
the tree closure assumes away.

The theory therefore makes a falsifiable prediction: hold P(k), m and theta
fixed, vary ONLY the overlap, and lambda_c must not move. This module builds
exactly that family.

Construction
------------
Layer 1 is a configuration-model hypergraph. Layer 2 reproduces a randomly
chosen fraction o of the layer-1 groups VERBATIM, and wires the remaining
member-slots at random:

    node v sits in d_S(v) of the copied groups, so it still needs k - d_S(v)
    layer-2 slots, which are filled by a configuration model on the residual
    stubs.

Because every node's layer-1 degree is k, d_S(v) <= k automatically, so the
residual degrees are non-negative and every node ends with layer-2 degree
exactly k. The joint degree distribution P(k) and the group sizes m are thus
identical across the whole family -- o is the only thing that moves.

o = 0 is the independent configuration-model baseline (its measured pair
overlap is O(1/N), not exactly zero); o = 1 makes layer 2 a verbatim copy of
layer 1, so every physical group carries both layers' transmission at once.

Measured statistic
------------------
    o_pair = #{node pairs sharing a layer-1 group AND a layer-2 group}
             / #{node pairs sharing a layer-1 group}

reported alongside the construction parameter, since o_pair is what actually
counts the 4-cycles the tree closure discards.
"""
import numpy as np
from itertools import combinations


# --------------------------------------------------------------------------
# measurement
# --------------------------------------------------------------------------
def _pair_set(groups):
    s = set()
    for row in groups:
        r = sorted(int(x) for x in row)
        s.update(combinations(r, 2))
    return s


def overlap(groups1, groups2):
    """Fraction of layer-1 co-membership pairs that recur in layer 2."""
    p1 = _pair_set(groups1)
    if not p1:
        return 0.0
    p2 = _pair_set(groups2)
    return len(p1 & p2) / len(p1)


# --------------------------------------------------------------------------
# construction
# --------------------------------------------------------------------------
def _config_groups(stubs, msize, rng, tries=24):
    """Cut shuffled stubs into groups of size msize, avoiding a node appearing
    twice in one group (a plain shuffle-and-cut produces such groups at rate
    O(k/N); they would make the group size effectively smaller)."""
    stubs = stubs.copy()
    rng.shuffle(stubs)
    n = len(stubs) // msize
    g = [stubs[i * msize:(i + 1) * msize].copy() for i in range(n)]
    for _ in range(tries):
        bad = [i for i, r in enumerate(g) if len(set(int(x) for x in r)) < msize]
        if not bad:
            break
        # repair by swapping a duplicated slot with a random slot elsewhere
        for i in bad:
            r = g[i]
            seen = {}
            for pos in range(msize):
                v = int(r[pos])
                if v in seen:
                    j = int(rng.integers(len(g)))
                    if j == i:
                        continue
                    q = int(rng.integers(msize))
                    r[pos], g[j][q] = g[j][q], r[pos]
                else:
                    seen[v] = pos
    g = [r for r in g if len(set(int(x) for x in r)) == msize]
    return g


def build_overlap_multiplex(N, m, k, o_target, rng):
    """Two-layer multiplex; every node has layer-degree (k, k); group sizes m;
    a fraction o_target of the layer-1 groups is reproduced verbatim in layer 2.

    Returns (groups, o_pair_realised).
    """
    m1, m2 = m
    if m1 != m2:
        raise ValueError("verbatim copying requires m1 == m2")
    g1 = _config_groups(np.repeat(np.arange(N), k), m1, rng)

    G = len(g1)
    n_copy = int(round(o_target * G))
    pick = rng.permutation(G)[:n_copy]
    copied = [g1[i].copy() for i in pick]

    used = np.zeros(N, dtype=np.int64)
    for r in copied:
        for v in r:
            used[int(v)] += 1
    resid = k - used
    if (resid < 0).any():                       # cannot happen for k1 == k2
        raise RuntimeError("negative residual degree")
    rest = _config_groups(np.repeat(np.arange(N), resid), m2, rng)

    g2 = copied + rest
    return [g1, g2], overlap(g1, g2)
