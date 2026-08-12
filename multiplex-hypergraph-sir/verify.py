"""Verification suite: reproduce every analytic anchor stated in the outline."""
import numpy as np
from theory import (cascade_C, next_gen_matrix, spectral_radius, rho_of_lambda,
                    lambda_c_rho, Closure)

def ok(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name}  {detail}")
    return cond

allpass = True

# ---- 1. C(2,lambda,1) = lambda/(1+lambda) = T  ------------------------------
for lam in [0.13, 0.35, 1.0, 2.0]:
    c = cascade_C(2, lam, 1)
    T = lam / (1 + lam)
    allpass &= ok(f"C(2,{lam},1)=T", abs(c - T) < 1e-12, f"C={c:.10f} T={T:.10f}")

# ---- 2. C(m,1,1) exceeds (m-1)T by 11/22/31%  (m=3,4,5) ---------------------
T1 = 0.5
for m, want in [(3, 0.11), (4, 0.22), (5, 0.31)]:
    c = cascade_C(m, 1.0, 1)
    excess = c / ((m - 1) * T1) - 1
    allpass &= ok(f"C({m},1,1) excess ~{want:.0%}", abs(excess - want) < 0.005,
                  f"C={c:.5f}  excess={excess:.4f}")

# ---- 3. theta>=2 => C=0 (single seed cannot activate) -----------------------
for m in [3, 4, 5]:
    allpass &= ok(f"C({m},1,theta=2)=0", cascade_C(m, 1.0, 2) == 0.0)

# ---- 4. ANCHOR A: P={(2,2),(3,3)}, m=(3,3), lam=0.13 -> rho=1.013 -----------
P = {(2, 2): 0.5, (3, 3): 0.5}
N = next_gen_matrix(P, m=(3, 3), lam=0.13)
rho = spectral_radius(N)
allpass &= ok("Anchor A  N11=N22~0.386", abs(N[0, 0] - 0.386) < 5e-4,
              f"N11={N[0,0]:.4f}")
allpass &= ok("Anchor A  rho~1.013", abs(rho - 1.013) < 5e-4, f"rho={rho:.4f}")

# ---- 5. Pairwise-degenerate: 5-regular single layer -> lambda_c = 1/3 -------
P5 = {(5,): 1.0}
lc = lambda_c_rho(P5, m=(2,))
allpass &= ok("Degenerate 5-regular lam_c=1/3", abs(lc - 1/3) < 1e-9,
              f"lam_c={lc:.10f}")

# ---- 6. ANCHOR B: theta>=2 channel does not move the threshold -------------
# Outline sec. 2.9: with P={(2,2),(3,3)}, m=(3,4), adding a theta=2 layer leaves
# the threshold essentially unchanged (a high-threshold channel changes outbreak
# *size*, not outbreak *onset*).  Under rho(N)=1 the theta=2 layer has C=0, so it
# drops out EXACTLY and lambda_c equals the single-(layer-1) value.
P1 = {(2,): 0.5, (3,): 0.5}                      # layer-1 marginal, m=3
lc_single = lambda_c_rho(P1, m=(3,))
lc_theta12 = lambda_c_rho(P, m=(3, 4), theta=(1, 2))
allpass &= ok("Anchor B  theta=2 layer drops out exactly",
              abs(lc_theta12 - lc_single) < 1e-9,
              f"lc(theta=1,2)={lc_theta12:.6f} == lc(single)={lc_single:.6f}")

# The outline quotes 0.410619 for this config -- that is the finite-t_max ODE
# *bisection* threshold, which sec. 3.1 documents as biased +2.5-2.6% high by the
# integration window.  rho(N)=1 gives the unbiased value; check the ratio.
bias = 0.410619 / lc_single - 1.0
allpass &= ok("Anchor B  outline 0.410619 = rho-value x (1 + 2.5%) [t_max bias]",
              0.024 < bias < 0.027, f"lc(rho)={lc_single:.6f}  bias=+{100*bias:.2f}%")

print("\n=== ALL PASS ===" if allpass else "\n!!! SOME FAILURES !!!")
