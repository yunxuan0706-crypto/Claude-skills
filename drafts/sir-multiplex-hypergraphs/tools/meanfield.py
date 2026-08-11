#!/usr/bin/env python3
"""Degree-based mean field for the multiplex hypergraph SIR, integrated in time.

The group activation probability is taken as the ensemble average
Theta_a = 1 - (1 - phi_a)^(m_a - 1), which is what discards the intra-group
correlation, the neighbourhood correlation and the within-group cascade all at
once. This is the baseline the closure is measured against.
"""
import math


class MeanField:
    def __init__(self, P, ms, betas):
        self.P, self.ms, self.betas, self.M = P, ms, betas, len(ms)
        self.ks = list(P)
        self.kbar = [sum(p*k[a] for k, p in P.items()) for a in range(self.M)]

    def run(self, eps, tmax, dt):
        n = len(self.ks)
        y = [1.0-eps]*n + [eps]*n            # s_k then i_k
        ms, betas, M, ks, P = self.ms, self.betas, self.M, self.ks, self.P

        def rhs(y):
            s, i = y[:n], y[n:]
            phi = [sum(P[k]*k[a]*i[j] for j, k in enumerate(ks))/self.kbar[a]
                   for a in range(M)]
            Th = [1.0-(1.0-phi[a])**(ms[a]-1) for a in range(M)]
            o = [0.0]*(2*n)
            for j, k in enumerate(ks):
                haz = sum(k[a]*betas[a]*Th[a] for a in range(M))
                o[j] = -haz*s[j]
                o[n+j] = haz*s[j] - i[j]      # mu = 1
            return o

        out, t = [], 0.0
        for _ in range(int(round(tmax/dt))+1):
            S = sum(P[k]*y[j] for j, k in enumerate(ks))
            I = sum(P[k]*y[n+j] for j, k in enumerate(ks))
            out.append((t, S, I))
            k1 = rhs(y)
            k2 = rhs([y[q]+.5*dt*k1[q] for q in range(2*n)])
            k3 = rhs([y[q]+.5*dt*k2[q] for q in range(2*n)])
            k4 = rhs([y[q]+dt*k3[q] for q in range(2*n)])
            y = [y[q]+dt/6*(k1[q]+2*k2[q]+2*k3[q]+k4[q]) for q in range(2*n)]
            t += dt
        return out
