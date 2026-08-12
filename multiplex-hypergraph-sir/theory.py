"""
Group-closure theory for SIR on multiplex hypergraphs.

Implements, from the research outline:
  - within-group cascade C(m, lambda, theta)          [Eq. (2.8)-(2.9)]
  - next-generation matrix N_ab                        [Eq. (2.10)]
  - epidemic threshold via rho(N) = 1                  [Prop. 3, Eq. (2.11)]
  - the group-closure ODE system                       [Eq. (2.2),(2.4)-(2.7)]
  - final epidemic size R(inf) from the closure

Conventions (mu-normalised time, so recovery rate = 1, lambda_a = beta_a/mu):
  A degree distribution P is a dict  {k_vector(tuple): probability}.
  k_vector[a] = number of layer-a groups a node belongs to.
"""
from functools import lru_cache
import numpy as np
from itertools import product as iproduct


# --------------------------------------------------------------------------
# 1. Within-group cascade  C(m, lambda, theta)      Eq. (2.8)-(2.9)
# --------------------------------------------------------------------------
def cascade_C(m, lam, theta=1):
    """Expected number of *additional* infections a single seed causes inside
    an isolated m-member group with activation threshold theta.

    u(i,s): expected further infections starting from i infected, s susceptible.
        u(i,s) = a/(a+i) [1 + u(i+1,s-1)] + i/(a+i) u(i-1,s),  a = lam*s*1[i>=theta]
        u(0,s) = u(i,0) = 0
    C = u(1, m-1).
    Well-founded: both transitions strictly decrease 2s+i.
    """
    @lru_cache(maxsize=None)
    def u(i, s):
        if i == 0 or s == 0:
            return 0.0
        a = lam * s if i >= theta else 0.0
        denom = a + i
        return (a / denom) * (1.0 + u(i + 1, s - 1)) + (i / denom) * u(i - 1, s)

    val = u(1, m - 1)
    u.cache_clear()
    return val


# --------------------------------------------------------------------------
# 2. Joint-degree moments
# --------------------------------------------------------------------------
def degree_moments(P, M):
    """Return (mean_k, cross) with mean_k[a]=<k_a>, cross[a,b]=<k_a k_b>."""
    mean_k = np.zeros(M)
    cross = np.zeros((M, M))
    tot = 0.0
    for k, p in P.items():
        tot += p
        k = np.asarray(k, dtype=float)
        mean_k += p * k
        cross += p * np.outer(k, k)
    if not np.isclose(tot, 1.0):
        raise ValueError(f"degree distribution not normalised (sum={tot})")
    return mean_k, cross


# --------------------------------------------------------------------------
# 3. Next-generation matrix and threshold      Eq. (2.10), Prop. 3
# --------------------------------------------------------------------------
def next_gen_matrix(P, m, lam, w=None, theta=None):
    """N_ab = C(m_b, lam_b, theta_b) * ( <k_a k_b>/<k_a> - delta_ab ),
    with lam_b = lam * w_b."""
    M = len(m)
    if w is None:
        w = np.ones(M)
    if theta is None:
        theta = np.ones(M, dtype=int)
    mean_k, cross = degree_moments(P, M)
    Cb = np.array([cascade_C(m[b], lam * w[b], theta[b]) for b in range(M)])
    N = np.empty((M, M))
    for a in range(M):
        for b in range(M):
            excess = cross[a, b] / mean_k[a] - (1.0 if a == b else 0.0)
            N[a, b] = Cb[b] * excess
    return N


def spectral_radius(N):
    return float(np.max(np.abs(np.linalg.eigvals(N))))


def rho_of_lambda(lam, P, m, w=None, theta=None):
    return spectral_radius(next_gen_matrix(P, m, lam, w, theta))


def lambda_c_rho(P, m, w=None, theta=None, lo=1e-4, hi=50.0, tol=1e-12):
    """Threshold lambda_c solving rho(N(lambda)) = 1 by bisection.
    rho is monotone increasing in lambda under the common scaling lambda_a=lam*w_a."""
    f = lambda L: rho_of_lambda(L, P, m, w, theta) - 1.0
    flo, fhi = f(lo), f(hi)
    if flo > 0:
        raise ValueError("rho already > 1 at lo; lower lo.")
    if fhi < 0:
        raise ValueError("rho still < 1 at hi; raise hi.")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if fm < 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


# --------------------------------------------------------------------------
# 4. Generating functions on the support of P
# --------------------------------------------------------------------------
class GenFun:
    """Psi and psi_a (and gradients) for a discrete joint-degree distribution."""

    def __init__(self, P, M):
        self.M = M
        self.ks = np.array([k for k in P.keys()], dtype=float)      # (K, M)
        self.p = np.array([P[k] for k in P.keys()], dtype=float)    # (K,)
        self.mean_k, _ = degree_moments(P, M)

    def Psi(self, x):
        # sum_k P(k) prod_c x_c^{k_c}
        return float(np.sum(self.p * np.prod(x ** self.ks, axis=1)))

    def psi(self, a, x):
        # (1/<k_a>) sum_k k_a P(k) prod_c x_c^{k_c - delta_ac}
        ex = self.ks.copy()
        ex[:, a] -= 1.0
        return float(np.sum(self.ks[:, a] * self.p * np.prod(x ** ex, axis=1))
                     / self.mean_k[a])

    def dpsi(self, a, c, x):
        # d psi_a / d x_c
        ex = self.ks.copy()
        ex[:, a] -= 1.0
        expo = ex[:, c]
        base = self.ks[:, a] * self.p
        # derivative wrt x_c: multiply by expo * x_c^{expo-1}, keep other factors
        xsafe = x.copy()
        terms = np.ones(len(self.p))
        for cc in range(self.M):
            if cc == c:
                # expo * x_c^{expo-1}
                with np.errstate(divide="ignore", invalid="ignore"):
                    fac = np.where(expo == 0, 0.0, expo * xsafe[cc] ** (expo - 1))
                terms *= fac
            else:
                terms *= xsafe[cc] ** ex[:, cc]
        return float(np.sum(base * terms) / self.mean_k[a])


# --------------------------------------------------------------------------
# 5. Group-closure ODE system            Eq. (2.4)-(2.7), (2.2)
# --------------------------------------------------------------------------
class Closure:
    """The dimension-Sum_a C(m_a+1,2)+... group-closure ODE.

    State layout: concatenation over layers of x^{(a)}_{sir} on the simplex
    s+i+r = m_a-1.  R(inf) is obtained directly from the fixed point:
        R(inf) = 1 - (1-eps) Psi(Phi(inf)).
    """

    def __init__(self, P, m, lam, w=None, theta=None, eps=1e-3):
        self.M = len(m)
        self.m = list(m)
        self.eps = eps
        if w is None:
            w = np.ones(self.M)
        if theta is None:
            theta = np.ones(self.M, dtype=int)
        self.beta = np.array([lam * w[a] for a in range(self.M)])  # mu=1
        self.theta = list(theta)
        self.gf = GenFun(P, self.M)

        # enumerate (s,i,r) states per layer with s+i+r = m_a-1
        self.states = []
        self.index = []
        off = 0
        self.offsets = []
        for a in range(self.M):
            n = self.m[a] - 1
            st = [(s, i, r) for s in range(n + 1) for i in range(n + 1 - s)
                  for r in [n - s - i]]
            idx = {sir: j for j, sir in enumerate(st)}
            self.states.append(st)
            self.index.append(idx)
            self.offsets.append(off)
            off += len(st)
        self.dim = off

    # ---- packing helpers ----
    def initial(self):
        y = np.zeros(self.dim)
        for a in range(self.M):
            n = self.m[a] - 1
            for j, (s, i, r) in enumerate(self.states[a]):
                if r == 0:
                    from math import comb
                    y[self.offsets[a] + j] = (comb(n, i)
                                              * (1 - self.eps) ** s
                                              * self.eps ** i)
        return y

    def phi(self, y):
        Phi = np.empty(self.M)
        PhiA = np.empty(self.M)
        for a in range(self.M):
            xa = y[self.offsets[a]: self.offsets[a] + len(self.states[a])]
            Phi[a] = xa.sum()
            PhiA[a] = sum(xa[j] for j, (s, i, r) in enumerate(self.states[a])
                          if i >= self.theta[a])
        return Phi, PhiA

    def hazard(self, Phi, PhiA):
        h = np.empty(self.M)
        for a in range(self.M):
            psi_a = self.gf.psi(a, Phi)
            acc = 0.0
            for c in range(self.M):
                acc += self.beta[c] * PhiA[c] * self.gf.dpsi(a, c, Phi) / psi_a
            h[a] = acc
        return h

    def rhs(self, t, y):
        dy = np.zeros_like(y)
        Phi, PhiA = self.phi(y)
        h = self.hazard(Phi, PhiA)
        for a in range(self.M):
            ba, tha = self.beta[a], self.theta[a]
            off = self.offsets[a]
            idx = self.index[a]
            for j, (s, i, r) in enumerate(self.states[a]):
                x = y[off + j]
                # loss: S->I from (s,i,r) at rate (h+ba*act(i))*s -> (s-1,i+1,r)
                act_i = 1.0 if i >= tha else 0.0
                rate_inf_out = (h[a] + ba * act_i) * s
                # loss: recovery I->R at rate i -> (s,i-1,r+1)   (mu=1)
                rate_rec_out = i
                # loss: transmit to u at rate ba*act(i) (absorbing)
                rate_trans = ba * act_i
                dy[off + j] -= (rate_inf_out + rate_rec_out + rate_trans) * x
                # gains into (s,i,r):
                # from (s+1,i-1,r) via S->I with source-i = i-1
                if (s + 1, i - 1, r) in idx and i - 1 >= 0:
                    src = y[off + idx[(s + 1, i - 1, r)]]
                    act_src = 1.0 if (i - 1) >= tha else 0.0
                    dy[off + j] += (h[a] + ba * act_src) * (s + 1) * src
                # from (s,i+1,r-1) via recovery
                if (s, i + 1, r - 1) in idx and r - 1 >= 0:
                    src = y[off + idx[(s, i + 1, r - 1)]]
                    dy[off + j] += (i + 1) * src
        return dy

    # ---- observables ----
    def Rinf(self, tmax=2000.0, rtol=1e-10, atol=1e-12):
        from scipy.integrate import solve_ivp
        y0 = self.initial()
        sol = solve_ivp(self.rhs, [0, tmax], y0, method="LSODA",
                        rtol=rtol, atol=atol, dense_output=False)
        yend = sol.y[:, -1]
        Phi, _ = self.phi(yend)
        S = (1 - self.eps) * self.gf.Psi(Phi)
        return 1.0 - S, Phi, sol

    def jacobian_dfe(self):
        """Numerical Jacobian of the ODE at the disease-free fixed point
        (eps=0, all mass in all-susceptible state)."""
        import numpy as np
        y0 = np.zeros(self.dim)
        for a in range(self.M):
            # all-susceptible compartment (m_a-1, 0, 0)
            j = self.index[a][(self.m[a] - 1, 0, 0)]
            y0[self.offsets[a] + j] = 1.0
        n = self.dim
        J = np.zeros((n, n))
        d = 1e-7
        f0 = self.rhs(0.0, y0)
        for kdim in range(n):
            yp = y0.copy()
            yp[kdim] += d
            J[:, kdim] = (self.rhs(0.0, yp) - f0) / d
        return J
