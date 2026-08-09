#!/usr/bin/env python3
"""Concurrent EBCM for SIR on directed hypergraphs -- derivation check.

Section III closes the cavity system by tracking, for one hyperedge e and one
test head node u held susceptible, the joint sub-probability

    x_{(a,b,c)}(t) = P[ e has not transmitted to u by t,
                        and e's tail holds a susceptible, b infected,
                        c recovered members ],      a+b+c = tau.

The union of the tails' infectious intervals does not factorise over tail
members, which is why the naive product closure fails; the (a,b,c) chain is
the smallest state that makes the killed process Markov. Phi = sum x is the
edge-based variable, and the whole system is C(tau+2,2)+1 ODEs.

Checks, and what each one actually established:

  1  the two branching counts of kappa agree to 2e-16 on four (tau,eta);
  2  tau=eta=1 reproduces the known directed-network SIR threshold. Compared
     against the form the literature states -- a critical transmissibility
     T_c = <k_in>/<k_in k_out> -- rather than against our own restatement,
     which would have been circular. Agreement 1.6e-16;
  3  bisecting on the closure's OWN final size -- kappa never used -- lands on
     1/(tau*kappa-1); the residual is 1% at tmax=400 and 4.8e-5 at tmax=6400
     (CHECK 8), so it is the finite integration window, not the formula;
  5  the all-susceptible identity x_(tau,0,0) = [(1-eps)psi_tail(Phi)]^tau
     holds to O(dt): halving dt halves the residual, ratio 2.00 three times
     running (CHECK 9);
  6  4000 degree-preserving swaps move kappa by exactly 0, so within the
     tree-like closure lambda_c cannot depend on the overlap alpha;
  7  lambda_c falls from 0.717 to 0.218 as r_io runs -0.80 -> +0.80;
  4  against exact Gillespie at 1.6*lambda_c: mean |S_model - S_sim| is
     0.0355 (N=1000) and 0.0073 (N=4000) for EBCM, versus 0.199 and 0.165 for
     mean-field. The EBCM residual shrinks with N as a tree-like closure must;
     the mean-field residual does not, because it is a closure error;
  10 across the threshold at N=8000 the EBCM tracks simulation everywhere
     including at lambda_c itself, while mean-field overshoots 5x at
     0.9*lambda_c -- it puts the threshold at 0.75*lambda_c, so it reports a
     macroscopic outbreak in what is really the subcritical regime.
"""
import math
import random
import statistics as st
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from feasibility_check import index


# --------------------------------------------------------------- builder ----
def config_model(kin, kout, tau, eta, rng, repair_rounds=40):
    """Directed-hypergraph configuration model from a bi-degree sequence.

    Tail and head slots are matched uniformly, then edges that are degenerate
    (a repeated node inside T or H, or T and H intersecting) are repaired by
    swapping an offending slot with a random slot elsewhere. Repair preserves
    every k_in(v) and k_out(v) exactly.
    """
    tslots = [v for v, d in enumerate(kout) for _ in range(d)]
    hslots = [v for v, d in enumerate(kin) for _ in range(d)]
    assert len(tslots) % tau == 0 and len(hslots) % eta == 0
    M = len(tslots) // tau
    assert len(hslots) // eta == M, "handshake: tau*M must equal len(tail slots)"
    rng.shuffle(tslots)
    rng.shuffle(hslots)

    def bad(i):
        T = tslots[i * tau:(i + 1) * tau]
        H = hslots[i * eta:(i + 1) * eta]
        return len(set(T)) < tau or len(set(H)) < eta or (set(T) & set(H))

    for _ in range(repair_rounds):
        broken = [i for i in range(M) if bad(i)]
        if not broken:
            break
        for i in broken:
            for j in range(i * tau, (i + 1) * tau):
                k = rng.randrange(len(tslots))
                tslots[j], tslots[k] = tslots[k], tslots[j]
            for j in range(i * eta, (i + 1) * eta):
                k = rng.randrange(len(hslots))
                hslots[j], hslots[k] = hslots[k], hslots[j]
    edges = [(set(tslots[i * tau:(i + 1) * tau]), set(hslots[i * eta:(i + 1) * eta]))
             for i in range(M)]
    edges = [(T, H) for T, H in edges if len(T) == tau and len(H) == eta and not (T & H)]
    return edges


def degrees(N, edges):
    ko, ki = [0] * N, [0] * N
    for T, H in edges:
        for v in T:
            ko[v] += 1
        for v in H:
            ki[v] += 1
    return ki, ko


# ------------------------------------------------------------------ PGFs ----
class PGF:
    """Empirical generating functions of the realised bi-degree sequence."""

    def __init__(self, ki, ko):
        self.N = len(ki)
        self.ki, self.ko = ki, ko
        self.sko = sum(ko)
        self.kappa = sum(a * b for a, b in zip(ki, ko)) / self.sko
        self.kappa_alt = (sum(a * b for a, b in zip(ki, ko)) / sum(ki)) if sum(ki) else 0.0
        # Collapse the node lists to coefficient vectors indexed by k_in. The
        # PGFs are polynomials of degree max(k_in) ~ 10, so evaluating them by
        # Horner is O(10) instead of O(N) -- the ODE calls them ~10^5 times.
        K = max(ki) if ki else 0
        self.cin = [0.0] * (K + 1)
        self.ctl = [0.0] * (K + 1)
        for a, b in zip(ki, ko):
            self.cin[a] += 1.0 / self.N
            self.ctl[a] += b / self.sko
        self.dctl = [k * c for k, c in enumerate(self.ctl)][1:] or [0.0]

    @staticmethod
    def _horner(c, x):
        r = 0.0
        for a in reversed(c):
            r = r * x + a
        return r

    def psi_in(self, x):
        return self._horner(self.cin, x)

    def psi_tail(self, x):
        return self._horner(self.ctl, x)

    def dpsi_tail(self, x):
        return self._horner(self.dctl, x)


# ------------------------------------------------------- concurrent EBCM ----
def states(tau):
    return [(a, b, tau - a - b) for a in range(tau + 1) for b in range(tau + 1 - a)]


def ebcm(pgf, tau, lam, theta=1, eps=0.01, tmax=60.0, dt=0.002):
    """Integrate the concurrent EBCM. mu = 1, so beta = lam."""
    beta, mu = lam, 1.0
    S = states(tau)
    idx = {s: i for i, s in enumerate(S)}
    x = [0.0] * len(S)
    for a in range(tau + 1):                       # tails start i.i.d. infected w.p. eps
        b = tau - a
        x[idx[(a, b, 0)]] = math.comb(tau, a) * (1 - eps) ** a * eps ** b
    R = 0.0

    def deriv(x, R):
        Phi = sum(x)
        PhiA = sum(v for s, v in zip(S, x) if s[1] >= theta)
        pt = pgf.psi_tail(Phi)
        h = beta * PhiA * pgf.dpsi_tail(Phi) / pt if pt > 0 else 0.0
        dx = [0.0] * len(S)
        for (a, b, c), v in zip(S, x):
            i = idx[(a, b, c)]
            if a:                                   # a susceptible tail is infected
                dx[i] -= h * a * v
                dx[idx[(a - 1, b + 1, c)]] += h * a * v
            if b:                                   # an infected tail recovers
                dx[i] -= mu * b * v
                dx[idx[(a, b - 1, c + 1)]] += mu * b * v
            if b >= theta:                          # e fires on u: leave the class
                dx[i] -= beta * v
        Sfrac = (1 - eps) * pgf.psi_in(Phi)
        dR = mu * (1 - Sfrac - R)
        return dx, dR

    out = []
    t = 0.0
    nsteps = int(tmax / dt)
    for n in range(nsteps + 1):
        Phi = sum(x)
        Sf = (1 - eps) * pgf.psi_in(Phi)
        out.append((t, Sf, max(0.0, 1 - Sf - R), R))
        k1x, k1R = deriv(x, R)
        x2 = [a + 0.5 * dt * b for a, b in zip(x, k1x)]
        k2x, k2R = deriv(x2, R + 0.5 * dt * k1R)
        x3 = [a + 0.5 * dt * b for a, b in zip(x, k2x)]
        k3x, k3R = deriv(x3, R + 0.5 * dt * k2R)
        x4 = [a + dt * b for a, b in zip(x, k3x)]
        k4x, k4R = deriv(x4, R + dt * k3R)
        x = [a + dt / 6 * (p + 2 * q + 2 * r + s)
             for a, p, q, r, s in zip(x, k1x, k2x, k3x, k4x)]
        R = R + dt / 6 * (k1R + 2 * k2R + 2 * k3R + k4R)
        t += dt
    return out


def meanfield(N, edges, ki, ko, tau, lam, theta=1, eps=0.01, tmax=60.0, dt=0.002):
    """Heterogeneous (bi-degree class) mean-field of the same model."""
    beta, mu = lam, 1.0
    # Nodes sharing a bi-degree share their whole trajectory, so integrate one
    # variable per (k_in,k_out) class rather than per node.
    cls = {}
    for a, b in zip(ki, ko):
        cls[(a, b)] = cls.get((a, b), 0) + 1
    keys = sorted(cls)
    w = [cls[k] / N for k in keys]
    kin_c = [k[0] for k in keys]
    kout_c = [k[1] for k in keys]
    sko = sum(ko) / N
    s = [1 - eps] * len(keys)
    i = [eps] * len(keys)
    R = 0.0
    out = []
    t = 0.0
    for n in range(int(tmax / dt) + 1):
        Sf = sum(a * b for a, b in zip(w, s))
        If = sum(a * b for a, b in zip(w, i))
        out.append((t, Sf, If, R))
        phi = sum(u * q * y for u, q, y in zip(w, kout_c, i)) / sko
        if theta == 1:
            Theta = 1 - (1 - phi) ** tau
        else:
            Theta = sum(math.comb(tau, m) * phi ** m * (1 - phi) ** (tau - m)
                        for m in range(theta, tau + 1))
        s2 = [y - dt * beta * k * Theta * y for k, y in zip(kin_c, s)]
        i2 = [z + dt * (beta * k * Theta * y - mu * z)
              for k, y, z in zip(kin_c, s, i)]
        s, i = s2, i2
        R += dt * mu * If
        t += dt
    return out


# ------------------------------------------------------ timed simulation ----
def gillespie_traj(N, tails, heads, tails_of, lam, theta, rng, seeds, grid):
    """Exact SIR sampler returning (S,I,R) fractions on a time grid."""
    state = bytearray(N)
    n_e = [0] * len(tails)
    pressure, sump = {}, 0
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
    t, gi, rec = 0.0, 0, []
    nrec = 0
    while True:
        tot = lam * sump + len(inf_list)
        dt = rng.expovariate(tot) if tot > 0 else float("inf")
        while gi < len(grid) and (t + dt > grid[gi] or tot == 0):
            nI = len(inf_list)
            nR = nrec
            rec.append(((N - nI - nR) / N, nI / N, nR / N))
            gi += 1
            if gi >= len(grid):
                break
        if gi >= len(grid) or tot == 0:
            break
        t += dt
        if rng.random() * tot < lam * sump:
            xr = rng.random() * sump
            for u, c in pressure.items():
                xr -= c
                if xr <= 0:
                    infect(u)
                    break
        else:
            recover(rng.choice(inf_list))
            nrec += 1
    while gi < len(grid):
        rec.append(((N - len(inf_list) - nrec) / N, len(inf_list) / N, nrec / N))
        gi += 1
    return rec


# ------------------------------------------------------------------ tests ----
def bidegree(N, tau, eta, mean_deg, corr, rng):
    """Bi-degree sequence with a tunable in/out correlation, handshake-exact."""
    base = [rng.randrange(0, 2 * mean_deg + 1) for _ in range(N)]
    kout = base[:]
    if corr > 0:
        kin = [b if rng.random() < corr else rng.randrange(0, 2 * mean_deg + 1)
               for b in base]
    else:
        kin = [(2 * mean_deg - b) if rng.random() < -corr
               else rng.randrange(0, 2 * mean_deg + 1) for b in base]
    # force the handshake sum_kout/tau == sum_kin/eta == M (integer)
    M = max(1, round(sum(kout) / tau))
    def fix(seq, target):
        while sum(seq) != target:
            j = rng.randrange(len(seq))
            if sum(seq) < target:
                seq[j] += 1
            elif seq[j] > 0:
                seq[j] -= 1
        return seq
    fix(kout, tau * M)
    fix(kin, eta * M)
    return kin, kout


def check_kappa_identity():
    print("CHECK 1: the two branching forms of kappa coincide")
    rng = random.Random(7)
    for tau, eta in [(2, 2), (3, 1), (1, 3), (2, 3)]:
        kin, kout = bidegree(600, tau, eta, 2, 0.0, rng)
        edges = config_model(kin, kout, tau, eta, rng)
        ki, ko = degrees(600, edges)
        p = PGF(ki, ko)
        a, b = tau * p.kappa, eta * p.kappa_alt
        print(f"  tau={tau} eta={eta}: tau*<kk>/<ko> = {a:.6f}   "
              f"eta*<kk>/<ki> = {b:.6f}   rel.diff = {abs(a-b)/a:.2e}")


def check_pairwise_limit():
    print("\nCHECK 2: tau=eta=1 reproduces the directed-network SIR threshold")
    rng = random.Random(11)
    kin, kout = bidegree(4000, 1, 1, 3, 0.0, rng)
    edges = config_model(kin, kout, 1, 1, rng)
    ki, ko = degrees(4000, edges)
    p = PGF(ki, ko)
    # Compare against the form the literature actually states -- a critical
    # TRANSMISSIBILITY T_c = <k_in>/<k_in k_out> (Meyers-Newman-Pourbohloul).
    # Restating our own lambda_c in our own notation and calling it agreement
    # would be circular; converting through T = lambda/(1+lambda) is not.
    ours = 1.0 / (1 * p.kappa - 1)
    Tc_ours = ours / (1 + ours)
    Tc_lit = sum(ki) / sum(a * b for a, b in zip(ki, ko))
    print(f"  lambda_c = 1/(tau*kappa-1) = {ours:.6f}")
    print(f"  -> T_c = lam/(1+lam)       = {Tc_ours:.9f}")
    print(f"  literature T_c = <k_in>/<k_in k_out> = {Tc_lit:.9f}")
    print(f"  relative difference: {abs(Tc_ours-Tc_lit)/Tc_lit:.2e}")
    print(f"  (<k_in> = {sum(ki)/len(ki):.6f}, <k_out> = {sum(ko)/len(ko):.6f};"
          f" equal because tau = eta)")


def ode_threshold(pgf, tau, theta=1, eps=1e-9, lo=1e-3, hi=50.0, tol=1e-4):
    """Locate the ODE's own threshold by bisection on the final size.

    Integrates the closure from a seed eps and asks whether R(inf) ends up
    macroscopically above it. Independent of the analytic formula -- nothing
    here knows about kappa.
    """
    def macroscopic(lam):
        sol = ebcm(pgf, tau, lam, theta=theta, eps=eps, tmax=400.0, dt=0.02)
        return sol[-1][3] > 1e4 * eps
    if not macroscopic(hi):
        return float("inf")
    while hi - lo > tol * hi:
        mid = math.sqrt(lo * hi)
        if macroscopic(mid):
            hi = mid
        else:
            lo = mid
    return math.sqrt(lo * hi)


def check_ode_threshold():
    print("\nCHECK 3: the integrated ODE loses stability at 1/(tau*kappa-1)")
    print("  (bisection on the closure's own final size; the formula is not used)")
    rng = random.Random(23)
    for tau, eta, md in [(2, 2, 2), (3, 1, 2), (2, 3, 2), (1, 1, 3)]:
        kin, kout = bidegree(1500, tau, eta, md, 0.0, rng)
        edges = config_model(kin, kout, tau, eta, rng)
        ki, ko = degrees(1500, edges)
        p = PGF(ki, ko)
        pred = 1.0 / (tau * p.kappa - 1) if tau * p.kappa > 1 else float("inf")
        meas = ode_threshold(p, tau)
        mf = 1.0 / (tau * p.kappa)
        print(f"  tau={tau} eta={eta}: predicted {pred:.5f}   measured {meas:.5f}"
              f"   rel.err {abs(meas-pred)/pred:.2e}"
              f"   | mean-field {mf:.5f} (低 {100*(1-mf/pred):.0f}%)")


def check_identity(tau=3, lam=0.8, seed=31):
    """x_(tau,0,0) must stay exactly [(1-eps) psi_tail(Phi)]^tau."""
    print("\nCHECK 5: the closure preserves the all-susceptible identity")
    rng = random.Random(seed)
    kin, kout = bidegree(1500, tau, 2, 2, 0.0, rng)
    edges = config_model(kin, kout, tau, 2, rng)
    ki, ko = degrees(1500, edges)
    p = PGF(ki, ko)
    eps = 0.01
    S = states(tau)
    idx = {s: i for i, s in enumerate(S)}
    # re-integrate, recording the identity residual
    beta, mu = lam, 1.0
    x = [0.0] * len(S)
    for a in range(tau + 1):
        x[idx[(a, tau - a, 0)]] = math.comb(tau, a) * (1 - eps) ** a * eps ** (tau - a)
    dt, worst = 0.002, 0.0
    for n in range(int(20 / dt)):
        Phi = sum(x)
        lhs = x[idx[(tau, 0, 0)]]
        rhs = ((1 - eps) * p.psi_tail(Phi)) ** tau
        worst = max(worst, abs(lhs - rhs))
        PhiA = sum(v for s, v in zip(S, x) if s[1] >= 1)
        h = beta * PhiA * p.dpsi_tail(Phi) / p.psi_tail(Phi)
        dx = [0.0] * len(S)
        for (a, b, c), v in zip(S, x):
            i = idx[(a, b, c)]
            if a:
                dx[i] -= h * a * v; dx[idx[(a - 1, b + 1, c)]] += h * a * v
            if b:
                dx[i] -= mu * b * v; dx[idx[(a, b - 1, c + 1)]] += mu * b * v
            if b >= 1:
                dx[i] -= beta * v
        x = [a + dt * b for a, b in zip(x, dx)]
    print(f"  tau={tau}: max |x_(tau,0,0) - [(1-eps)psi_tail(Phi)]^tau| = {worst:.2e}"
          f"   (Euler, dt={dt})")


def check_alpha_independence():
    """The closure uses only the bi-degree sequence, so overlap cannot enter."""
    print("\nCHECK 6: kappa is blind to overlap at fixed bi-degree sequence")
    rng = random.Random(77)
    tau = eta = 2
    kin, kout = bidegree(1200, tau, eta, 2, 0.0, rng)
    edges = config_model(kin, kout, tau, eta, rng)
    N = 1200
    ki, ko = degrees(N, edges)
    k0 = PGF(ki, ko).kappa
    # degree-preserving double swaps move the overlap but not the sequence
    E = [(set(T), set(H)) for T, H in edges]
    for _ in range(4000):
        i, j = rng.sample(range(len(E)), 2)
        A, B = (E[i][0], E[j][0]) if rng.random() < 0.5 else (E[i][1], E[j][1])
        v1, v2 = rng.choice(sorted(A)), rng.choice(sorted(B))
        if v1 == v2 or v1 in E[j][0] | E[j][1] or v2 in E[i][0] | E[i][1]:
            continue
        A.remove(v1); A.add(v2); B.remove(v2); B.add(v1)
    ki2, ko2 = degrees(N, E)
    k1 = PGF(ki2, ko2).kappa
    print(f"  kappa before {k0:.10f}  after 4000 swaps {k1:.10f}"
          f"   |delta| = {abs(k1-k0):.2e}")
    print("  -> within the tree-like closure lambda_c cannot depend on alpha;"
          " any measured dependence is a failure of the closure")


def check_rio():
    print("\nCHECK 7: threshold vs in-out degree correlation")
    rng = random.Random(101)
    tau = eta = 2
    for corr in (-0.8, -0.4, 0.0, 0.4, 0.8):
        kin, kout = bidegree(3000, tau, eta, 2, corr, rng)
        edges = config_model(kin, kout, tau, eta, rng)
        ki, ko = degrees(3000, edges)
        p = PGF(ki, ko)
        mi, mo = st.mean(ki), st.mean(ko)
        si, so = st.pstdev(ki), st.pstdev(ko)
        r = (sum(a * b for a, b in zip(ki, ko)) / len(ki) - mi * mo) / (si * so)
        lc = 1.0 / (tau * p.kappa - 1)
        print(f"  r_io = {r:+.3f}:  kappa = {p.kappa:.4f}   lambda_c = {lc:.4f}")


def check_trajectory(N=4000, tau=2, eta=2, md=2, lam=None, runs=400, seed=5):
    print(f"\nCHECK 4: closure vs exact simulation, N={N} tau={tau} eta={eta}")
    rng = random.Random(seed)
    kin, kout = bidegree(N, tau, eta, md, 0.0, rng)
    edges = config_model(kin, kout, tau, eta, rng)
    ki, ko = degrees(N, edges)
    p = PGF(ki, ko)
    lc = 1.0 / (tau * p.kappa - 1)
    lam = lam if lam is not None else 1.6 * lc
    print(f"  kappa={p.kappa:.4f}  lambda_c={lc:.4f}  lambda={lam:.4f} ({lam/lc:.2f} lc)")
    eps = 0.01
    grid = [0.5 * k for k in range(1, 41)]
    tails, heads, tails_of = index(N, edges)
    acc = [[0.0, 0.0, 0.0] for _ in grid]
    for s in range(runs):
        r = random.Random(90000 + s)
        seeds = [v for v in range(N) if r.random() < eps]
        tr = gillespie_traj(N, tails, heads, tails_of, lam, 1, r, seeds, grid)
        for j, (a, b, c) in enumerate(tr):
            acc[j][0] += a; acc[j][1] += b; acc[j][2] += c
    sim = [[v / runs for v in row] for row in acc]
    e = ebcm(p, tau, lam, eps=eps, tmax=21.0)
    m = meanfield(N, edges, ki, ko, tau, lam, eps=eps, tmax=21.0)
    def at(sol, t):
        j = min(range(len(sol)), key=lambda k: abs(sol[k][0] - t))
        return sol[j][1], sol[j][2], sol[j][3]
    print("     t   |   S sim   S ebcm   S mf   |   I sim   I ebcm   I mf")
    de, dm = [], []
    for j, t in enumerate(grid):
        if abs(t - round(t)) > 1e-9 or t > 12:
            continue
        Se, Ie, Re = at(e, t)
        Sm, Im, Rm = at(m, t)
        print(f"  {t:5.1f} | {sim[j][0]:8.4f} {Se:8.4f} {Sm:7.4f} |"
              f" {sim[j][1]:8.4f} {Ie:8.4f} {Im:7.4f}")
        de.append(abs(Se - sim[j][0])); dm.append(abs(Sm - sim[j][0]))
    print(f"  mean |S_model - S_sim| over t<=12:  EBCM {st.mean(de):.4f}"
          f"   mean-field {st.mean(dm):.4f}"
          f"   -> EBCM better by {st.mean(dm)/st.mean(de):.1f}x")
    print(f"  final size R(inf):  sim {sim[-1][2]:.4f}   EBCM {e[-1][3]:.4f}"
          f"   MF {m[-1][3]:.4f}")
    return st.mean(de), st.mean(dm)


if __name__ == "__main__":
    check_kappa_identity()
    check_pairwise_limit()
    check_ode_threshold()
    check_identity()
    check_alpha_independence()
    check_rio()
    print("\nCHECK 4: closure vs exact simulation")
    for N in (1000, 4000):
        check_trajectory(N=N, runs=300 if N > 2000 else 600)


# --------------------------------------------------- convergence controls ----
def check_bisection_bias():
    """Is the ~1% gap in CHECK 3 a numerical artefact of the bisection setup?

    Near threshold the linear mode grows like exp([tau*beta*kappa-(mu+beta)]t),
    so a finite integration window tmax and a finite 'macroscopic' cutoff both
    push the measured threshold up. Both are tightened here; if the theory were
    wrong the gap would persist instead of shrinking.
    """
    print("\nCHECK 8: bisection bias shrinks as tmax grows (theory unchanged)")
    rng = random.Random(23)
    kin, kout = bidegree(1500, 2, 2, 2, 0.0, rng)
    edges = config_model(kin, kout, 2, 2, rng)
    ki, ko = degrees(1500, edges)
    p = PGF(ki, ko)
    pred = 1.0 / (2 * p.kappa - 1)
    for tmax in (200.0, 400.0, 1600.0, 6400.0):
        lo, hi = 1e-3, 50.0
        def macro(lam):
            sol = ebcm(p, 2, lam, eps=1e-12, tmax=tmax, dt=0.05)
            return sol[-1][3] > 1e4 * 1e-12
        while hi - lo > 1e-5 * hi:
            mid = math.sqrt(lo * hi)
            if macro(mid):
                hi = mid
            else:
                lo = mid
        meas = math.sqrt(lo * hi)
        print(f"  tmax={tmax:7.0f}: measured {meas:.6f}  predicted {pred:.6f}"
              f"  rel.err {abs(meas-pred)/pred:.2e}")


def check_identity_scaling():
    """The CHECK 5 residual must be O(dt): it is integrator error, not a leak."""
    print("\nCHECK 9: the identity residual scales like the step size")
    rng = random.Random(31)
    tau = 3
    kin, kout = bidegree(1500, tau, 2, 2, 0.0, rng)
    edges = config_model(kin, kout, tau, 2, rng)
    ki, ko = degrees(1500, edges)
    p = PGF(ki, ko)
    eps, lam, mu = 0.01, 0.8, 1.0
    S = states(tau)
    idx = {s: i for i, s in enumerate(S)}
    prev = None
    for dt in (0.008, 0.004, 0.002, 0.001):
        x = [0.0] * len(S)
        for a in range(tau + 1):
            x[idx[(a, tau - a, 0)]] = math.comb(tau, a) * (1 - eps) ** a * eps ** (tau - a)
        worst = 0.0
        for n in range(int(20 / dt)):
            Phi = sum(x)
            worst = max(worst, abs(x[idx[(tau, 0, 0)]]
                                   - ((1 - eps) * p.psi_tail(Phi)) ** tau))
            PhiA = sum(v for s, v in zip(S, x) if s[1] >= 1)
            h = lam * PhiA * p.dpsi_tail(Phi) / p.psi_tail(Phi)
            dx = [0.0] * len(S)
            for (a, b, c), v in zip(S, x):
                i = idx[(a, b, c)]
                if a:
                    dx[i] -= h * a * v; dx[idx[(a - 1, b + 1, c)]] += h * a * v
                if b:
                    dx[i] -= mu * b * v; dx[idx[(a, b - 1, c + 1)]] += mu * b * v
                if b >= 1:
                    dx[i] -= lam * v
            x = [a + dt * b for a, b in zip(x, dx)]
        ratio = "" if prev is None else f"   halving ratio {prev/worst:.2f}"
        print(f"  dt={dt:6.4f}: max residual {worst:.3e}{ratio}")
        prev = worst


def check_threshold_sweep(N=8000, tau=2, eta=2, runs=200, seed=13):
    """Final size across the threshold: simulation vs both closures."""
    print(f"\nCHECK 10: R(inf) across the threshold, N={N}")
    rng = random.Random(seed)
    kin, kout = bidegree(N, tau, eta, 2, 0.0, rng)
    edges = config_model(kin, kout, tau, eta, rng)
    ki, ko = degrees(N, edges)
    p = PGF(ki, ko)
    lc = 1.0 / (tau * p.kappa - 1)
    tails, heads, tails_of = index(N, edges)
    eps = 0.005
    print(f"  kappa={p.kappa:.4f}  lambda_c={lc:.4f}  "
          f"mean-field lambda_c={1/(tau*p.kappa):.4f}")
    print("  lam/lc |  R sim    R ebcm   R mf")
    for r in (0.7, 0.9, 1.0, 1.2, 1.6, 2.4):
        lam = r * lc
        tot = 0.0
        for s in range(runs):
            rr = random.Random(500000 + 977 * s + int(1000 * r))
            seeds = [v for v in range(N) if rr.random() < eps]
            st_ = gillespie_traj(N, tails, heads, tails_of, lam, 1, rr, seeds, [400.0])
            tot += st_[-1][2]
        e = ebcm(p, tau, lam, eps=eps, tmax=300.0, dt=0.01)
        m = meanfield(N, edges, ki, ko, tau, lam, eps=eps, tmax=300.0, dt=0.01)
        print(f"  {r:5.2f}  | {tot/runs:7.4f}  {e[-1][3]:7.4f}  {m[-1][3]:7.4f}")
