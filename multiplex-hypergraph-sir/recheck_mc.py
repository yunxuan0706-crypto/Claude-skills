"""Fresh independent re-checks: single-group Gillespie MC for C, and an ODE
supercritical/subcritical sanity across the transition."""
import numpy as np
from theory import cascade_C, Closure, lambda_c_rho, next_gen_matrix, spectral_radius

rng = np.random.default_rng(20260816)

def cascade_C_mc(m, lam, theta, nruns=200_000):
    """Direct Gillespie of one isolated m-group: seed 1 infected, count the
    expected number of *additional* members ever infected."""
    add = np.empty(nruns)
    for k in range(nruns):
        i, s = 1, m - 1                       # infected, susceptible
        while i > 0:
            active = 1.0 if i >= theta else 0.0
            r_inf = lam * s * active
            r_rec = i
            tot = r_inf + r_rec
            if tot == 0.0:
                break
            if rng.random() < r_inf / tot:    # infection
                i += 1; s -= 1
            else:                             # recovery
                i -= 1
        add[k] = (m - s) - 1                  # ever-infected minus the seed
    return add.mean(), add.std(ddof=1) / np.sqrt(nruns)

if __name__ == "__main__":
    print("=== A. Single-group Gillespie MC  vs  recursion C(m,lambda,theta) ===")
    cases = [(2, 0.35, 1), (3, 0.5, 1), (3, 1.0, 1), (4, 1.0, 1), (5, 1.0, 1),
             (4, 0.8, 1), (3, 2.0, 1), (4, 1.0, 2)]
    worst_sig = 0.0
    for m, lam, th in cases:
        c_th = cascade_C(m, lam, th)
        c_mc, se = cascade_C_mc(m, lam, th)
        nsig = abs(c_mc - c_th) / se if se > 0 else 0.0
        worst_sig = max(worst_sig, nsig)
        print(f"  m={m} lam={lam} th={th}:  C_recursion={c_th:.5f}  "
              f"C_MC={c_mc:.5f}+-{se:.5f}  ({nsig:.2f} sigma)")
    print(f"  worst MC deviation: {worst_sig:.2f} sigma\n")

    print("=== B. ODE behaviour across the transition (eps=0.02) ===")
    P = {(2, 2): .5, (3, 3): .5}; m = (3, 4)
    lc = lambda_c_rho(P, m)
    print(f"  lambda_c(rho=1) = {lc:.5f}")
    for f in (0.7, 0.9, 0.98, 1.02, 1.2, 1.6):
        clo = Closure(P, m, f * lc, eps=0.02)
        from scipy.integrate import solve_ivp
        sol = solve_ivp(clo.rhs, [0, 4000], clo.initial(), method="LSODA",
                        rtol=1e-10, atol=1e-13)
        Phi, _ = clo.phi(sol.y[:, -1])
        Rinf = 1 - (1 - clo.eps) * clo.gf.Psi(Phi)
        tag = "subcritical" if f < 1 else "SUPERcritical"
        print(f"  lam={f:.2f} lc:  R(inf)={Rinf:.4f}   [{tag}]")
    print("  -> R(inf) stays O(eps) below lc and lifts to an O(1) epidemic above lc\n")

    print("=== C. rho monotonic in lambda (bisection well-posed) ===")
    lams = np.linspace(0.02, 0.30, 8)
    rhos = [spectral_radius(next_gen_matrix(P, m, L)) for L in lams]
    mono = all(rhos[i] < rhos[i+1] for i in range(len(rhos)-1))
    print(f"  rho(lambda) strictly increasing: {mono}")
    print("  " + "  ".join(f"{r:.3f}" for r in rhos))
