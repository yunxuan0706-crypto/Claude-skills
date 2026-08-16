"""
Degree mean-field closure and the term-by-term decomposition of its deviation.

Outline sec. 2.8 (Prop. 5). The mean-field next-generation matrix (2.13)

    N^MF_ab = (m_b - 1) * lam_b * <k^a k^b> / <k^a>

differs from the exact group-level matrix (2.10)

    N_ab = C(m_b, lam_b, th_b) * [ <k^a k^b>/<k^a> - delta_ab ]

in exactly three places, with two competing directions:

  (T) transmissibility : lam_b   instead of  T_b = lam_b/(1+lam_b).
      Mean field treats the whole infectious period as beta/mu, ignoring that a
      group is "used up" on a member once it has transmitted to it.
      Factor  f_T = T/lam = 1/(1+lam) < 1  ->  MF OVERestimates R.

  (D) excess subtraction : <k^a k^b>/<k^a>  instead of  ... - delta_ab.
      Mean field lets the infection turn back along the group it arrived by.
      Factor  f_D = rho(X - delta)/rho(X) < 1  ->  MF OVERestimates R.

  (C) cascade factor : (m_b-1) T_b  instead of  C(m_b, lam_b, th_b).
      Linearising member-by-member cannot see that intra-group infections
      lengthen the group's active period.
      Factor  f_C = C / [(m-1) T] > 1  ->  MF UNDERestimates R.

Written multiplicatively, each switch s in {0,1}:

    N(s_T, s_D, s_C)_ab = (m_b-1) lam_b * f_T(b)^{s_T} * f_C(b)^{s_C}
                          * [ X_ab - s_D * delta_ab ]

so (0,0,0) is the mean field and (1,1,1) is exact.
"""
import numpy as np

from theory import cascade_C, degree_moments, spectral_radius


# --------------------------------------------------------------------------
# switchable next-generation matrix
# --------------------------------------------------------------------------
def ngm_switched(P, m, lam, s_T=1, s_D=1, s_C=1, w=None, theta=None):
    """Next-generation matrix with each mean-field deviation switched on/off.

    s_T=s_D=s_C=1 -> exact (2.10);  s_T=s_D=s_C=0 -> mean field (2.13).
    """
    M = len(m)
    if w is None:
        w = np.ones(M)
    if theta is None:
        theta = np.ones(M, dtype=int)
    mean_k, cross = degree_moments(P, M)

    pref = np.empty(M)
    for b in range(M):
        lb = lam * w[b]
        Tb = lb / (1.0 + lb)
        base = (m[b] - 1) * lb                       # mean-field prefactor
        f_T = (Tb / lb) if lb > 0 else 1.0           # = 1/(1+lam_b)  < 1
        denom = (m[b] - 1) * Tb
        Cb = cascade_C(m[b], lb, theta[b])
        f_C = (Cb / denom) if denom > 0 else 1.0     # >= 1  (=0 when theta>=2)
        pref[b] = base * (f_T ** s_T) * (f_C ** s_C)

    N = np.empty((M, M))
    for a in range(M):
        for b in range(M):
            X = cross[a, b] / mean_k[a]
            N[a, b] = pref[b] * (X - (s_D if a == b else 0.0))
    return N


def rho_switched(P, m, lam, **kw):
    return spectral_radius(ngm_switched(P, m, lam, **kw))


def lambda_c_switched(P, m, s_T=1, s_D=1, s_C=1, w=None, theta=None,
                      lo=1e-6, hi=200.0, tol=1e-13):
    """Threshold for a given switch setting, by bisection on rho = 1."""
    f = lambda L: rho_switched(P, m, L, s_T=s_T, s_D=s_D, s_C=s_C,
                               w=w, theta=theta) - 1.0
    if f(lo) > 0 or f(hi) < 0:
        raise ValueError("threshold not bracketed")
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


def lambda_c_mf(P, m, **kw):
    return lambda_c_switched(P, m, s_T=0, s_D=0, s_C=0, **kw)


# --------------------------------------------------------------------------
# the three factors, reported individually
# --------------------------------------------------------------------------
def factors(P, m, lam, w=None, theta=None):
    """Return (f_T, f_D, f_C, net) as they act on the spectral radius.

    net = rho_exact / rho_MF.  net < 1 means mean field overestimates R.
    f_D is evaluated spectrally (it is not a per-layer scalar in general).
    """
    M = len(m)
    if w is None:
        w = np.ones(M)
    if theta is None:
        theta = np.ones(M, dtype=int)
    lam_b = lam * np.asarray(w, float)
    T_b = lam_b / (1.0 + lam_b)
    f_T = float(np.mean(1.0 / (1.0 + lam_b)))                    # per-layer, = T/lam
    C_b = np.array([cascade_C(m[b], lam_b[b], theta[b]) for b in range(M)])
    with np.errstate(divide="ignore", invalid="ignore"):
        fc = np.where((np.asarray(m) - 1) * T_b > 0,
                      C_b / ((np.asarray(m) - 1) * T_b), 1.0)
    f_C = float(np.mean(fc))
    # spectral effect of the excess subtraction alone
    r_with = rho_switched(P, m, lam, s_T=0, s_D=1, s_C=0, w=w, theta=theta)
    r_without = rho_switched(P, m, lam, s_T=0, s_D=0, s_C=0, w=w, theta=theta)
    f_D = r_with / r_without
    net = (rho_switched(P, m, lam, s_T=1, s_D=1, s_C=1, w=w, theta=theta)
           / r_without)
    return f_T, f_D, f_C, net


# --------------------------------------------------------------------------
# mean-field ODE (2.12) and its disease-free Jacobian
# --------------------------------------------------------------------------
class MeanFieldODE:
    """Degree mean field:  sdot_k = -s_k sum_a k_a beta_a Theta_a,
                           idot_k =  s_k sum_a k_a beta_a Theta_a - i_k,
       with Theta_a = 1 - (1 - phi_a)^{m_a - 1},
            phi_a   = sum_k k_a P(k) i_k / <k^a>          (mu = 1)."""

    def __init__(self, P, m, lam, w=None, theta=None):
        self.M = len(m)
        self.m = list(m)
        if w is None:
            w = np.ones(self.M)
        self.beta = lam * np.asarray(w, float)
        self.ks = np.array(list(P.keys()), dtype=float)     # (K, M)
        self.p = np.array([P[k] for k in P.keys()], dtype=float)
        self.mean_k, _ = degree_moments(P, self.M)
        self.K = len(self.p)

    def phi(self, i_k):
        return (self.ks * self.p[:, None] * i_k[:, None]).sum(0) / self.mean_k

    def rhs(self, t, y):
        s_k = y[:self.K]
        i_k = y[self.K:]
        ph = self.phi(i_k)
        Theta = 1.0 - (1.0 - ph) ** (np.asarray(self.m) - 1)
        force = (self.ks * self.beta * Theta).sum(1)         # per degree class
        ds = -s_k * force
        di = s_k * force - i_k
        return np.concatenate([ds, di])

    def jacobian_dfe(self, d=1e-6):
        """Numerical Jacobian of the i_k block at the disease-free state.

        Central differences: the forward difference carries an O(d) truncation
        error that grows with m (Theta = 1-(1-phi)^{m-1} has curvature
        ~ (m-1)(m-2)), which would swamp the O(1e-7) agreement being tested.
        Note phi is perturbed to +-d, so the state leaves the physical
        simplex on the minus side; that is harmless here because the rhs is a
        smooth polynomial in phi.
        """
        y0 = np.concatenate([np.ones(self.K), np.zeros(self.K)])
        J = np.zeros((self.K, self.K))
        for j in range(self.K):
            yp, ym = y0.copy(), y0.copy()
            yp[self.K + j] += d
            ym[self.K + j] -= d
            J[:, j] = (self.rhs(0.0, yp) - self.rhs(0.0, ym))[self.K:] / (2 * d)
        return J

    def spectral_abscissa(self):
        return float(np.max(np.real(np.linalg.eigvals(self.jacobian_dfe()))))
