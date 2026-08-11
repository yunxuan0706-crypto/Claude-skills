#!/usr/bin/env python3
"""Does the group-count closure actually reproduce the stochastic process?

The threshold check (tools/verify_threshold.py) only tests the linearisation.
This tests the full nonlinear equations against exact Gillespie simulation on
multiplex hypergraphs, and measures how the error scales with N -- which is what
decides whether Prop 3.1 ("exact in the tree limit") may be claimed at all.
"""
import math, random, sys
from collections import defaultdict

# ------------------------------------------------------------ generator ----
def config_model(kseq, m, rng, repair=200):
    """kseq[v] = number of groups of size m containing v."""
    slots = [v for v, k in enumerate(kseq) for _ in range(k)]
    rng.shuffle(slots)
    n = len(slots) // m
    groups = [slots[i*m:(i+1)*m] for i in range(n)]
    for _ in range(repair):
        bad = [gi for gi, g in enumerate(groups) if len(set(g)) < m]
        if not bad:
            break
        for gi in bad:
            for pos in range(m):
                gj = rng.randrange(len(groups)); pj = rng.randrange(m)
                groups[gi][pos], groups[gj][pj] = groups[gj][pj], groups[gi][pos]
    return [g for g in groups if len(set(g)) == m]

def build(P, ms, N, rng):
    ks = list(P)
    counts = [max(1, int(round(P[k]*N))) for k in ks]
    kvec = []
    for k, c in zip(ks, counts):
        kvec += [k]*c
    rng.shuffle(kvec)
    N = len(kvec)
    layers = []
    for a, m in enumerate(ms):
        layers.append(config_model([kv[a] for kv in kvec], m, rng))
    return N, layers

# ------------------------------------------------------------ Gillespie ----
def gillespie(N, layers, betas, tmax, dt, eps, rng):
    """Exact CTMC sampling; mu = 1. theta = 1 for every layer."""
    S, I, R = 0, 1, 2
    state = [S]*N
    ginfo, of = [], [[] for _ in range(N)]          # groups, and each node's groups
    for a, gs in enumerate(layers):
        for g in gs:
            gid = len(ginfo)
            ginfo.append([a, g, len(g), 0])          # layer, members, #S, #I
            for v in g:
                of[v].append(gid)
    buckets = defaultdict(list)                      # (a,s) -> [gid]
    pos = {}
    def bkey(gid):
        a, g, s, i = ginfo[gid]
        return (a, s) if (i >= 1 and s >= 1) else None
    def brefresh(gid, old):
        new = bkey(gid)
        if old == new: return
        if old is not None:
            lst = buckets[old]; p = pos[gid]; last = lst.pop()
            if p < len(lst): lst[p] = last; pos[last] = p
            del pos[gid]
        if new is not None:
            buckets[new].append(gid); pos[gid] = len(buckets[new])-1
    inf = []
    for v in range(N):
        if rng.random() < eps:
            state[v] = I; inf.append(v)
    for v in inf:
        for gid in of[v]:
            old = bkey(gid); ginfo[gid][2] -= 1; ginfo[gid][3] += 1; brefresh(gid, old)
    t, out, nxt = 0.0, [], 0.0
    while t < tmax:
        while nxt <= t and nxt <= tmax:
            out.append((nxt, state.count(S)/N)); nxt += dt
        W = sum(betas[a]*s*len(l) for (a, s), l in buckets.items() if l)
        Rr = float(len(inf))
        tot = W + Rr
        if tot <= 0: break
        t += rng.expovariate(tot)
        if rng.random() < W/tot:                      # infection
            x = rng.random()*W; pick = None
            for (a, s), l in buckets.items():
                if not l: continue
                w = betas[a]*s*len(l)
                if x < w: pick = l[int(x/w*len(l)) % len(l)]; break
                x -= w
            if pick is None: continue
            g = ginfo[pick][1]
            cand = [v for v in g if state[v] == S]
            if not cand: continue
            v = rng.choice(cand)
            state[v] = I; inf.append(v)
            for gid in of[v]:
                old = bkey(gid); ginfo[gid][2] -= 1; ginfo[gid][3] += 1; brefresh(gid, old)
        else:                                          # recovery
            j = rng.randrange(len(inf)); v = inf[j]
            inf[j] = inf[-1]; inf.pop()
            state[v] = R
            for gid in of[v]:
                old = bkey(gid); ginfo[gid][3] -= 1; brefresh(gid, old)
    while nxt <= tmax:
        out.append((nxt, state.count(S)/N)); nxt += dt
    return out

# ---------------------------------------------------------------- closure --
class Closure:
    def __init__(self, P, ms, betas):
        self.P, self.ms, self.betas, self.M = P, ms, betas, len(ms)
        self.kbar = [sum(p*k[a] for k, p in P.items()) for a in range(self.M)]
        self.st = [[(s, i, m-1-s-i) for s in range(m) for i in range(m-s)] for m in ms]
        self.ix = [{v: j for j, v in enumerate(S)} for S in self.st]
        self.off = [0]
        for S in self.st: self.off.append(self.off[-1]+len(S))
    def Psi(self, Phi):
        return sum(p*math.prod(Phi[c]**k[c] for c in range(self.M)) for k, p in self.P.items())
    def psi_a(self, a, Phi):
        tot, d = 0.0, [0.0]*self.M
        for k, p in self.P.items():
            if k[a] == 0: continue
            ex = [k[c]-(1 if c == a else 0) for c in range(self.M)]
            w = p*k[a]
            tot += w*math.prod(Phi[c]**ex[c] for c in range(self.M))
            for c in range(self.M):
                if ex[c] > 0:
                    d[c] += w*ex[c]*Phi[c]**(ex[c]-1)*math.prod(
                        Phi[e]**ex[e] for e in range(self.M) if e != c)
        return tot/self.kbar[a], [x/self.kbar[a] for x in d]
    def rhs(self, y):
        M = self.M; n = self.off[-1]
        x, Rv = y[:n], y[n]
        Phi = [sum(x[self.off[a]:self.off[a+1]]) for a in range(M)]
        PhiA = [sum(x[self.off[a]+j] for j, (s, i, r) in enumerate(self.st[a]) if i >= 1)
                for a in range(M)]
        h = []
        for a in range(M):
            ps, dp = self.psi_a(a, Phi)
            h.append(sum(self.betas[c]*PhiA[c]*dp[c]/ps for c in range(M)))
        o = [0.0]*(n+1)
        for a in range(M):
            b, ix = self.betas[a], self.ix[a]
            for (s, i, r), j in ix.items():
                v = x[self.off[a]+j]
                if v == 0.0: continue
                rate = (h[a] + (b if i >= 1 else 0.0))*s
                if rate:
                    o[self.off[a]+j] -= rate*v
                    o[self.off[a]+ix[(s-1, i+1, r)]] += rate*v
                if i:
                    o[self.off[a]+j] -= i*v
                    o[self.off[a]+ix[(s, i-1, r+1)]] += i*v
                    o[self.off[a]+j] -= b*v
        Sv = (1-self.eps)*self.Psi(Phi)
        o[n] = 1.0*(1.0 - Sv - Rv)                    # dR/dt = mu I
        return o
    def run(self, eps, tmax, dt):
        self.eps = eps
        n = self.off[-1]
        y = [0.0]*(n+1)
        for a in range(self.M):
            m = self.ms[a]
            for i in range(m):
                s = m-1-i
                y[self.off[a]+self.ix[a][(s, i, 0)]] = math.comb(m-1, i)*(1-eps)**s*eps**i
        out, t = [], 0.0
        nsteps = int(round(tmax/dt))
        for k in range(nsteps+1):
            Phi = [sum(y[self.off[a]:self.off[a+1]]) for a in range(self.M)]
            out.append((t, (1-eps)*self.Psi(Phi)))
            k1 = self.rhs(y)
            k2 = self.rhs([y[i]+.5*dt*k1[i] for i in range(len(y))])
            k3 = self.rhs([y[i]+.5*dt*k2[i] for i in range(len(y))])
            k4 = self.rhs([y[i]+dt*k3[i] for i in range(len(y))])
            y = [y[i]+dt/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]) for i in range(len(y))]
            t += dt
        return out, y

    def invariants(self, eps, tmax, dt):
        """sum rule  d(sum x)/dt = -beta PhiA, and the all-susceptible identity"""
        self.eps = eps
        _, y = self.run(eps, tmax, dt)
        n = self.off[-1]
        Phi = [sum(y[self.off[a]:self.off[a+1]]) for a in range(self.M)]
        worst_id = 0.0
        for a in range(self.M):
            lhs = y[self.off[a]+self.ix[a][(self.ms[a]-1, 0, 0)]]
            ps, _ = self.psi_a(a, Phi)
            rhs = ((1-eps)*ps)**(self.ms[a]-1)
            worst_id = max(worst_id, abs(lhs-rhs))
        f = self.rhs(y)
        worst_sum = 0.0
        for a in range(self.M):
            dsum = sum(f[self.off[a]:self.off[a+1]])
            PhiA = sum(y[self.off[a]+j] for j, (s, i, r) in enumerate(self.st[a]) if i >= 1)
            worst_sum = max(worst_sum, abs(dsum + self.betas[a]*PhiA))
        return worst_id, worst_sum


if __name__ == "__main__":
    P = {(2, 2): .5, (3, 3): .5}; ms = (3, 4); lam = 0.35
    betas = [lam, lam]; eps, tmax, dt = 0.02, 12.0, 0.005
    cl = Closure(P, ms, betas)
    ode, _ = cl.run(eps, tmax, dt)
    ode = {round(t, 6): s for t, s in ode}
    print(f"closure vs exact simulation   lambda={lam}, m=({ms[0]},{ms[1]}), eps={eps}")
    print(f"{'N':>7} {'runs':>5} {'mean |dS| over t<=12':>22} {'max |dS|':>10}")
    for N0, runs in ((500, 400), (1000, 300), (2000, 200), (4000, 120)):
        rng = random.Random(7 + N0)
        acc, cnt = None, 0
        for r in range(runs):
            N, layers = build(P, ms, N0, random.Random(1000 + N0*13 + r))
            tr = gillespie(N, layers, betas, tmax, 0.25, eps, rng)
            if acc is None: acc = [0.0]*len(tr)
            for j, (_, s) in enumerate(tr): acc[j] += s
            cnt += 1
        sim = [(j*0.25, acc[j]/cnt) for j in range(len(acc))]
        d = [abs(s - ode[round(t, 6)]) for t, s in sim if t <= tmax]
        print(f"{N0:>7} {runs:>5} {sum(d)/len(d):22.5f} {max(d):10.5f}")
