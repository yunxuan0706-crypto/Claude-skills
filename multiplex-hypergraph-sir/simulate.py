"""
Exact stochastic simulation of SIR on a multiplex hypergraph.

This is the ground truth for Figure 3: the mean-field / group-closure
comparison there is otherwise theory-against-theory, so the direction and size
of the deviation rests on the closure being right. Here the same quantity is
measured by direct continuous-time (Gillespie) simulation of the microscopic
process on a configuration-model multiplex hypergraph -- with real loops and
real finite-N effects, no closure of any kind.

Observable: the mean total outbreak size from a single random seed in the
subcritical regime, which branching theory predicts to be

    chi(lambda) = 1^T (I - N)^{-1} v_0 .

Comparing the simulated chi against the exact N and against the mean-field
N^MF tests the three-factor decomposition directly, including its sign.
"""
import numpy as np


# --------------------------------------------------------------------------
# configuration-model multiplex hypergraph
# --------------------------------------------------------------------------
def build_multiplex(P, m, N, rng):
    """Return groups[a] = list of node-index arrays, and the realised degrees.

    Each node draws a layer-degree vector from P; layer-a stubs are shuffled
    and cut into groups of size m_a (the remainder is discarded, which is an
    O(m/N) effect on the degree sequence).
    """
    M = len(m)
    keys = list(P.keys())
    probs = np.array([P[k] for k in keys], float)
    idx = rng.choice(len(keys), size=N, p=probs)
    kmat = np.array([keys[i] for i in idx], dtype=np.int64)      # (N, M)

    groups = []
    for a in range(M):
        stubs = np.repeat(np.arange(N), kmat[:, a])
        rng.shuffle(stubs)
        n_grp = len(stubs) // m[a]
        g = stubs[:n_grp * m[a]].reshape(n_grp, m[a])
        groups.append(g)
    return groups, kmat


def _membership(groups, N):
    """node -> list of (layer, group index) it belongs to, as flat CSR arrays."""
    M = len(groups)
    owner = [[] for _ in range(N)]
    for a in range(M):
        for gi, members in enumerate(groups[a]):
            for v in members:
                owner[v].append((a, gi))
    return owner


# --------------------------------------------------------------------------
# Gillespie SIR with group activation
# --------------------------------------------------------------------------
def outbreak_size(groups, owner, beta, theta, seed_node, rng, cap=None):
    """One outbreak from a single seed; returns the total number ever infected.

    Rates (mu = 1): every infected node recovers at rate 1; every group with at
    least theta_a infected members transmits to each of its susceptible members
    at rate beta_a. Only the infected set and the groups touching it are ever
    examined, so subcritical outbreaks cost O(outbreak size).
    """
    M = len(groups)
    INF, REC = 1, 2
    state = {}                       # node -> INF/REC  (absent = susceptible)
    state[seed_node] = INF
    infected = {seed_node}
    n_ever = 1
    # groups that currently contain >= 1 infected member
    touched = {}                     # (a, gi) -> number infected in it

    def touch(v, delta):
        for (a, gi) in owner[v]:
            touched[(a, gi)] = touched.get((a, gi), 0) + delta

    touch(seed_node, +1)

    while infected:
        # --- build the event list ---
        rec_rate = float(len(infected))
        trans = []                   # (rate, a, gi, susceptible members)
        tot_trans = 0.0
        for (a, gi), n_inf in touched.items():
            if n_inf < theta[a]:
                continue
            members = groups[a][gi]
            sus = [int(v) for v in members if v not in state]
            if not sus:
                continue
            r = beta[a] * len(sus)
            trans.append((r, a, gi, sus))
            tot_trans += r

        total = rec_rate + tot_trans
        if total <= 0:
            break
        u = rng.random() * total
        if u < rec_rate:                                  # a recovery
            # all infected recover at the same rate, but *which* one recovers
            # matters (they sit in different groups), so choose uniformly
            inf_list = list(infected)
            v = inf_list[rng.integers(len(inf_list))]
            infected.discard(v)
            state[v] = REC
            touch(v, -1)
        else:                                             # a transmission
            u -= rec_rate
            for (r, a, gi, sus) in trans:
                if u < r:
                    w = sus[rng.integers(len(sus))]
                    state[w] = INF
                    infected.add(w)
                    touch(w, +1)
                    n_ever += 1
                    break
                u -= r
        if cap is not None and n_ever >= cap:
            return n_ever, True                            # hit the cap
    return n_ever, False


def mean_outbreak_size(P, m, lam, N=20000, n_seeds=4000, n_graphs=3,
                       w=None, theta=None, seed=0, cap=None):
    """Simulated chi(lambda): mean total outbreak size from a single seed."""
    M = len(m)
    if w is None:
        w = np.ones(M)
    if theta is None:
        theta = [1] * M
    beta = [lam * w[a] for a in range(M)]
    if cap is None:
        cap = max(2000, N // 10)

    sizes, capped = [], 0
    rng = np.random.default_rng(seed)
    for g in range(n_graphs):
        groups, kmat = build_multiplex(P, m, N, rng)
        owner = _membership(groups, N)
        for _ in range(n_seeds // n_graphs):
            v = int(rng.integers(N))
            s, hit = outbreak_size(groups, owner, beta, theta, v, rng, cap=cap)
            sizes.append(s)
            capped += hit
    sizes = np.array(sizes, float)
    return dict(chi=sizes.mean(), sem=sizes.std(ddof=1) / np.sqrt(len(sizes)),
                n=len(sizes), capped=capped)
