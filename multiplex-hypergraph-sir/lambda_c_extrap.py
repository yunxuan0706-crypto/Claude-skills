"""
Determine lambda_c by extrapolation from subcritical final states.

Protocol (research outline, sec. 3.1):
  The mean subcritical outbreak size is  chi(lambda) = 1^T (I - N)^{-1} v0,
  which diverges as (lambda_c - lambda)^{-1} as rho(N) -> 1^-.  Hence the
  reciprocal  y(lambda) = eps / R(inf) = 1/chi  falls *linearly* to zero, and
  its zero-crossing is lambda_c.

  R(inf) here is the final epidemic size of the nonlinear group-closure ODE
  (Eq. 2.5), integrated to steady state -- an independent numerical path from
  the rho(N)=1 eigenvalue computation.

  Finite-seed bias is removed by Richardson extrapolation eps -> 0
  (R(inf)/eps is linear in eps for small eps).
"""
import numpy as np
from scipy.integrate import solve_ivp
from theory import Closure, lambda_c_rho


def _Rinf(clo, eps, tmax, rtol=1e-11, atol=1e-14):
    y0 = clo.initial()
    sol = solve_ivp(clo.rhs, [0, tmax], y0, method="LSODA", rtol=rtol, atol=atol)
    Phi, PhiA = clo.phi(sol.y[:, -1])
    S = (1 - eps) * clo.gf.Psi(Phi)
    leftover = float(PhiA.sum())          # activated-and-untransmitted mass at t=tmax
    return 1.0 - S, leftover


def chi_subcritical(P, m, lam, w=None, theta=None, eps=1e-4, tmax=60000.0):
    """eps->0 mean outbreak size chi = lim_{eps->0} R(inf)/eps, via Richardson."""
    clo1 = Closure(P, m, lam, w=w, theta=theta, eps=eps)
    clo2 = Closure(P, m, lam, w=w, theta=theta, eps=eps / 2)
    R1, lo1 = _Rinf(clo1, eps, tmax)
    R2, lo2 = _Rinf(clo2, eps / 2, tmax)
    leftover = max(lo1, lo2)
    f1, f2 = R1 / eps, R2 / (eps / 2)
    return (2 * f2 - f1), leftover


def _xintercept(x, y):
    """OLS y=a+bx; x0=-a/b with delta-method standard error."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    n = len(x)
    X = np.vstack([np.ones(n), x]).T
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    a, b = beta
    resid = y - X @ beta
    s2 = (resid @ resid) / (n - 2)
    cov = s2 * np.linalg.inv(X.T @ X)
    x0 = -a / b
    g = np.array([-1.0 / b, a / b**2])
    return x0, float(np.sqrt(g @ cov @ g)), (a, b)


def run_config(P, m, w=None, theta=None, eps=1e-4,
               frac=(0.90, 0.995), npts=10, tmax=60000.0):
    lc = lambda_c_rho(P, m, w=w, theta=theta)
    lams = lc * np.linspace(frac[0], frac[1], npts)
    ys, leftovers = [], []
    for lam in lams:
        chi, leftover = chi_subcritical(P, m, lam, w=w, theta=theta,
                                        eps=eps, tmax=tmax)
        ys.append(1.0 / chi)
        leftovers.append(leftover)
    ys = np.array(ys)
    x0, se, (a, b) = _xintercept(lams, ys)
    return dict(P=P, m=m, theta=theta, lc_rho=lc, lc_extrap=x0, sigma=se,
                dev=x0 - lc, dev_sigma=(x0 - lc) / se,
                lams=lams, ys=ys, fit=(a, b),
                max_leftover=float(np.max(leftovers)))


# canonical configuration set (spans lambda_c ~ 0.08 .. 0.40)
CONFIGS = [
    ("2-layer m=(4,4)",  {(2, 2): .5, (3, 3): .5}, (4, 4), None),
    ("2-layer m=(3,4)",  {(2, 2): .5, (3, 3): .5}, (3, 4), None),
    ("2-layer m=(3,3)",  {(2, 2): .5, (3, 3): .5}, (3, 3), None),
    ("3-reg. m=4",       {(3,): 1.0},              (4,),   None),
    ("2-layer m=(2,3)",  {(2, 2): .5, (3, 3): .5}, (2, 3), None),
    ("6-reg. pairwise",  {(6,): 1.0},              (2,),   None),
    ("5-reg. pairwise",  {(5,): 1.0},              (2,),   None),
    ("m=3 deg{2,3}",     {(2,): .5, (3,): .5},     (3,),   None),
]

if __name__ == "__main__":
    import json
    results = []
    print(f"{'config':18s} {'lc(rho=1)':>10s} {'lc(extrap)':>11s} "
          f"{'dev%':>8s} {'sigma':>9s} {'dev/sig':>8s} {'leftover':>9s}")
    for name, P, m, theta in CONFIGS:
        r = run_config(P, m, theta=theta)
        r["name"] = name
        results.append(r)
        print(f"{name:18s} {r['lc_rho']:10.6f} {r['lc_extrap']:11.6f} "
              f"{100*r['dev']/r['lc_rho']:+8.4f} {r['sigma']:9.2e} "
              f"{r['dev_sigma']:+8.2f} {r['max_leftover']:9.1e}")

    maxsig = max(abs(r["dev_sigma"]) for r in results)
    maxrel = max(abs(r["dev"] / r["lc_rho"]) for r in results)
    print(f"\nmax |dev| = {maxsig:.2f} sigma   max rel. dev. = {maxrel:.2e}")

    # persist for the figure
    out = []
    for r in results:
        out.append(dict(name=r["name"], lc_rho=r["lc_rho"],
                        lc_extrap=r["lc_extrap"], sigma=r["sigma"],
                        dev_sigma=r["dev_sigma"],
                        lams=r["lams"].tolist(), ys=r["ys"].tolist(),
                        fit=list(r["fit"])))
    with open("figure_data.json", "w") as f:
        json.dump(out, f, indent=2)
    print("wrote figure_data.json")
