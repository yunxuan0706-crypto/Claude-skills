"""Verify the group-closure ODE: sum rule (2.6), identity (2.7), DFE Jacobian."""
import numpy as np
from theory import Closure, lambda_c_rho, next_gen_matrix, spectral_radius

def ok(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name}  {detail}")
    return cond

allpass = True

# Six representative configurations (doc sec. 1.6): pairwise-degenerate,
# anti-correlated degrees, asymmetric rates and group sizes, three layers,
# and two standard multiplex cases. (theta>=2 drop-out is Anchor B in verify.py.)
configs = [
    ("pairwise 5-reg (lc=1/3)", {(5,): 1.0}, (2,), None, None),
    ("anti-correlated {2,4}",   {(2, 4): 0.5, (4, 2): 0.5}, (3, 4), None, None),
    ("asymmetric m=(3,5) w=(1,2)", {(2, 2): 0.5, (3, 3): 0.5}, (3, 5),
                                np.array([1.0, 2.0]), None),
    ("three-layer M=3",         {(2, 2, 2): 0.5, (3, 3, 3): 0.5}, (3, 3, 3), None, None),
    ("2-layer (3,4)",           {(2, 2): 0.5, (3, 3): 0.5}, (3, 4), None, None),
    ("single m=3",              {(2,): 0.5, (3,): 0.5}, (3,), None, None),
]

for name, P, m, w, theta in configs:
    M = len(m)
    lam = 0.30
    clo = Closure(P, m, lam, w=w, theta=theta, eps=0.02)
    y = clo.initial()

    # ---- sum rule (2.6): d/dt Phi_a == -beta_a Phi^A_a  ----
    # take a few points along an RK4-integrated trajectory
    from scipy.integrate import solve_ivp
    sol = solve_ivp(clo.rhs, [0, 8], y, method="LSODA", rtol=1e-11, atol=1e-13,
                    t_eval=[0.5, 1.0, 2.0, 4.0])
    maxres = 0.0
    for col in range(sol.y.shape[1]):
        yc = sol.y[:, col]
        dy = clo.rhs(0.0, yc)
        # dPhi_a from dy
        Phi, PhiA = clo.phi(yc)
        for a in range(M):
            off = clo.offsets[a]; nst = len(clo.states[a])
            dPhi_a = dy[off:off+nst].sum()
            res = abs(dPhi_a - (-clo.beta[a] * PhiA[a]))
            maxres = max(maxres, res)
    allpass &= ok(f"[{name}] sum rule (2.6)", maxres < 1e-10,
                  f"max residual={maxres:.2e}")

    # ---- identity (2.7): x_{m-1,0,0} = ((1-eps) psi_a(Phi))^{m-1}  ----
    maxres2 = 0.0
    for col in range(sol.y.shape[1]):
        yc = sol.y[:, col]
        Phi, _ = clo.phi(yc)
        for a in range(M):
            j = clo.index[a][(m[a]-1, 0, 0)]
            lhs = yc[clo.offsets[a] + j]
            rhs = ((1-clo.eps) * clo.gf.psi(a, Phi)) ** (m[a]-1)
            maxres2 = max(maxres2, abs(lhs - rhs))
    allpass &= ok(f"[{name}] identity (2.7)", maxres2 < 1e-8,
                  f"max residual={maxres2:.2e}")

    # ---- DFE Jacobian spectral abscissa vs rho(N)=1 threshold ----
    # The abscissa vanishes for EVERY lambda <= lambda_c (the recovered /
    # already-transmitted directions contribute structural zero eigenvalues),
    # so "abscissa == 0 at lambda_c" is one-sided: it cannot detect a threshold
    # that is too low. Bisect instead for the lambda at which the abscissa
    # first turns positive, and compare THAT with lambda_c.
    lc = lambda_c_rho(P, m, w=w, theta=theta)

    def abscissa(L):
        J = Closure(P, m, L, w=w, theta=theta, eps=0.0).jacobian_dfe()
        return float(np.max(np.real(np.linalg.eigvals(J))))

    lo, hi = 0.2 * lc, 3.0 * lc
    assert abscissa(lo) <= 1e-9 < abscissa(hi)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if abscissa(mid) > 1e-9:
            hi = mid
        else:
            lo = mid
    lc_jac = 0.5 * (lo + hi)
    # central-difference Jacobian (theory.py) matches rho(N)=1 to the bisection
    # resolution ~1e-9; the old forward difference floored at ~1e-7 and would
    # fail this bound -- a guard against regressing that.
    allpass &= ok(f"[{name}] Jacobian instability onset == lam_c",
                  abs(lc_jac / lc - 1) < 1e-8,
                  f"lam_c={lc:.8f}  onset={lc_jac:.8f}  rel={abs(lc_jac/lc-1):.1e}")
    allpass &= ok(f"[{name}]   and abscissa is 0 at lam_c, >0 just above",
                  abs(abscissa(lc)) < 1e-9 and abscissa(1.02 * lc) > 1e-6)

print("\n=== ALL PASS ===" if allpass else "\n!!! SOME FAILURES !!!")
