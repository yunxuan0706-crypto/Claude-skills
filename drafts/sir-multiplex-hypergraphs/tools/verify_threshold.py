"""Decisive check of Prop. 3.2: does linearising the group-count closure
reproduce rho(N)=1 with N_ab = C(m_b,lam_b)[<k^a k^b>/<k^a> - delta_ab]?

Two completely different routes: numerical Jacobian of the ODE system at the
disease-free state, versus a multitype branching count. They must agree.
"""
import numpy as np, itertools
from functools import lru_cache

def C(m, lam):
    @lru_cache(None)
    def u(i, s):
        if i == 0 or s == 0: return 0.0
        tot = lam*s + i
        return (lam*s/tot)*(1.0+u(i+1, s-1)) + (i/tot)*u(i-1, s)
    return u(1, m-1)

class Model:
    def __init__(self, P, ms, ws, lam):
        self.P = P                      # {(k1,k2): prob}
        self.ms = ms; self.ws = ws
        self.betas = [lam*w for w in ws]      # mu = 1
        self.M = len(ms)
        self.kbar = [sum(p*k[a] for k,p in P.items()) for a in range(self.M)]
        self.states = [[(s,i,m-1-s-i) for s in range(m) for i in range(m-s)]
                       for m in ms]
        self.idx = [{st:j for j,st in enumerate(S)} for S in self.states]
        self.off = np.cumsum([0]+[len(S) for S in self.states])

    def psi_a(self, a, Phi):
        """excess PGF for layer a, and its partial derivatives"""
        tot, d = 0.0, [0.0]*self.M
        for k, p in self.P.items():
            if k[a] == 0: continue
            w = p*k[a]
            base = 1.0
            expo = [k[c]-(1 if c==a else 0) for c in range(self.M)]
            for c in range(self.M): base *= Phi[c]**expo[c]
            tot += w*base
            for c in range(self.M):
                if expo[c] > 0:
                    d[c] += w*expo[c]*Phi[c]**(expo[c]-1)*np.prod(
                        [Phi[e]**expo[e] for e in range(self.M) if e != c])
        n = self.kbar[a]
        return tot/n, [x/n for x in d]

    def rhs(self, x):
        M, ms = self.M, self.ms
        Phi = [x[self.off[a]:self.off[a+1]].sum() for a in range(M)]
        PhiA = [sum(x[self.off[a]+j] for j,(s,i,r) in enumerate(self.states[a]) if i>=1)
                for a in range(M)]
        h = []
        for a in range(M):
            psi, dpsi = self.psi_a(a, Phi)
            h.append(sum(self.betas[c]*PhiA[c]*dpsi[c]/psi for c in range(M)))
        out = np.zeros_like(x)
        for a in range(M):
            b, ix = self.betas[a], self.idx[a]
            for (s,i,r), j in ix.items():
                v = x[self.off[a]+j]
                rate_out = (h[a] + b*(1 if i>=1 else 0))*s          # a susceptible gets infected
                out[self.off[a]+j] -= rate_out*v
                if s >= 1:
                    out[self.off[a]+ix[(s-1,i+1,r)]] += rate_out*v
                out[self.off[a]+j] -= 1.0*i*v                        # recovery (mu=1)
                if i >= 1:
                    out[self.off[a]+ix[(s,i-1,r+1)]] += 1.0*i*v
                if i >= 1:                                           # transmission to u (killing)
                    out[self.off[a]+j] -= b*v
        return out

    def dfe(self):
        x = np.zeros(self.off[-1])
        for a in range(self.M):
            x[self.off[a]+self.idx[a][(self.ms[a]-1,0,0)]] = 1.0
        return x

    def growth_rate(self):
        x0 = self.dfe(); n = len(x0); eps = 1e-7
        J = np.zeros((n,n))
        f0 = self.rhs(x0)
        for j in range(n):
            xp = x0.copy(); xp[j] += eps
            J[:,j] = (self.rhs(xp)-f0)/eps
        return max(np.linalg.eigvals(J).real)

    def N_matrix(self):
        Mn = np.zeros((self.M,self.M))
        for a in range(self.M):
            for b in range(self.M):
                kab = sum(p*k[a]*k[b] for k,p in self.P.items())
                Mn[a,b] = C(self.ms[b], self.betas[b])*(kab/self.kbar[a] - (1 if a==b else 0))
        return Mn

def lc_from(fn, lo=0.001, hi=5.0):
    for _ in range(80):
        mid = 0.5*(lo+hi)
        if fn(mid) > 0: hi = mid
        else: lo = mid
    return 0.5*(lo+hi)

cases = [
 ("uncorrelated, m=(3,4)", {(2,2):.25,(2,4):.25,(4,2):.25,(4,4):.25}, (3,4), (1.0,1.0)),
 ("positively corr.",      {(2,2):.5,(4,4):.5},                        (3,4), (1.0,1.0)),
 ("negatively corr.",      {(2,4):.5,(4,2):.5},                        (3,4), (1.0,1.0)),
 ("asymmetric weights",    {(2,2):.5,(4,4):.5},                        (3,5), (1.0,0.5)),
 ("pairwise m=(2,2)",      {(2,3):.5,(4,1):.5},                        (2,2), (1.0,1.0)),
]
print(f"{'case':26s} {'lc (ODE Jacobian)':>18s} {'lc (rho(N)=1)':>15s} {'rel.diff':>10s}")
for name, P, ms, ws in cases:
    f_ode = lambda L: Model(P, ms, ws, L).growth_rate()
    f_bp  = lambda L: max(abs(np.linalg.eigvals(Model(P,ms,ws,L).N_matrix())))-1.0
    a, b = lc_from(f_ode), lc_from(f_bp)
    print(f"{name:26s} {a:18.6f} {b:15.6f} {abs(a-b)/b:10.2e}")
