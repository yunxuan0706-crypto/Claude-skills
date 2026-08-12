"""Verify the group-closure ODE: sum rule (2.6), identity (2.7), DFE Jacobian."""
import numpy as np
from theory import Closure, lambda_c_rho, next_gen_matrix, spectral_radius

def ok(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name}  {detail}")
    return cond

allpass = True

# Configs to test
configs = [
    ("single m=3",     {(2,):0.5,(3,):0.5}, (3,), None, None),
    ("2-layer (3,4)",  {(2,2):0.5,(3,3):0.5}, (3,4), None, None),
    ("2-layer theta",  {(2,2):0.5,(3,3):0.5}, (3,4), None, (1,2)),
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
    # At lambda_c, the DFE Jacobian should have max real part eigenvalue = 0.
    lc = lambda_c_rho(P, m, w=w, theta=theta)
    clo_c = Closure(P, m, lc, w=w, theta=theta, eps=0.0)
    J = clo_c.jacobian_dfe()
    max_re = float(np.max(np.real(np.linalg.eigvals(J))))
    allpass &= ok(f"[{name}] DFE Jacobian abscissa=0 at lam_c", abs(max_re) < 1e-4,
                  f"lam_c={lc:.6f}  max Re eig={max_re:.2e}")

print("\n=== ALL PASS ===" if allpass else "\n!!! SOME FAILURES !!!")
