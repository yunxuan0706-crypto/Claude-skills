"""Verification of the mean-field closure and the three-term decomposition."""
import numpy as np
from theory import lambda_c_rho, cascade_C
from meanfield import (lambda_c_mf, lambda_c_switched, ngm_switched,
                       rho_switched, factors, MeanFieldODE)

def ok(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name}  {detail}")
    return cond

allpass = True
P = {(3, 3): 0.5, (5, 5): 0.5}          # k in {3,5}, two layers, m1=m2=m

# ---- 1. outline anchor: lam_c^MF / lam_c for m = 2,3,4,6,8 -----------------
want = {2: 0.76, 3: 0.85, 4: 0.88, 6: 0.90, 8: 0.91}
for m, w in want.items():
    r = lambda_c_mf(P, (m, m)) / lambda_c_rho(P, (m, m))
    allpass &= ok(f"lam_c^MF/lam_c at m={m} ~ {w}", abs(r - w) < 0.005,
                  f"ratio={r:.4f}")

# ---- 2. ratio is monotone increasing in m and never exceeds 1 -------------
ms = [2, 3, 4, 5, 6, 7, 8, 10, 12]
rs = [lambda_c_mf(P, (m, m)) / lambda_c_rho(P, (m, m)) for m in ms]
allpass &= ok("ratio < 1 for all m (MF always underestimates lam_c)",
              all(r < 1 for r in rs), f"max={max(rs):.4f}")
allpass &= ok("ratio monotone increasing in m",
              all(rs[i] < rs[i+1] for i in range(len(rs)-1)),
              "  ".join(f"{r:.3f}" for r in rs))

# ---- 3. switches: (1,1,1) == exact, (0,0,0) == mean field ------------------
for m in (2, 3, 5):
    ex = lambda_c_rho(P, (m, m))
    sw = lambda_c_switched(P, (m, m), s_T=1, s_D=1, s_C=1)
    allpass &= ok(f"switch(1,1,1) == exact at m={m}", abs(sw - ex) < 1e-9,
                  f"{sw:.10f} vs {ex:.10f}")
    mf = lambda_c_mf(P, (m, m))
    sw0 = lambda_c_switched(P, (m, m), s_T=0, s_D=0, s_C=0)
    allpass &= ok(f"switch(0,0,0) == mean field at m={m}", abs(sw0 - mf) < 1e-12)

# ---- 4. mean-field matrix vs mean-field ODE Jacobian -----------------------
# at lam_c^MF the disease-free Jacobian of (2.12) must have zero abscissa
for m in (2, 3, 4, 6):
    lc = lambda_c_mf(P, (m, m))
    ab = MeanFieldODE(P, (m, m), lc).spectral_abscissa()
    allpass &= ok(f"MF ODE Jacobian abscissa = 0 at lam_c^MF (m={m})",
                  abs(ab) < 1e-7, f"lam_c^MF={lc:.8f}  abscissa={ab:.2e}")

# ---- 5. factor directions: f_T<1, f_D<1, f_C>1, and they multiply ---------
for m in (3, 5):
    for lam in (0.05, 0.2, 1.0):
        fT, fD, fC, net = factors(P, (m, m), lam)
        prod = fT * fD * fC
        allpass &= ok(f"factors multiply to net (m={m}, lam={lam})",
                      abs(prod - net) < 1e-10,
                      f"f_T={fT:.4f} f_D={fD:.4f} f_C={fC:.4f} -> {prod:.6f} vs {net:.6f}")
        allpass &= ok(f"  directions  f_T<1<f_C, f_D<1 (m={m}, lam={lam})",
                      fT < 1 and fD < 1 and fC > 1)

# ---- 6. term-by-term: which correction dominates ---------------------------
print("\nTerm-by-term thresholds (P={(3,3),(5,5)}, m=(3,3)):")
m3 = (3, 3)
base = lambda_c_mf(P, m3)
ex = lambda_c_rho(P, m3)
names = {(0,0,0): "mean field (none on)",
         (1,0,0): "only T  (lam -> T)",
         (0,1,0): "only D  (excess -1)",
         (0,0,1): "only C  (cascade)",
         (1,1,0): "T + D",
         (1,0,1): "T + C",
         (0,1,1): "D + C",
         (1,1,1): "all on = exact"}
for s, nm in names.items():
    lc = lambda_c_switched(P, m3, s_T=s[0], s_D=s[1], s_C=s[2])
    print(f"  {nm:24s} lam_c={lc:.6f}   (x{lc/base:.4f} vs MF)")
print(f"  exact/MF ratio = {ex/base:.4f}")

# ---- 7. the net effect DOES change sign, but never at the threshold -------
# Outline sec. 2.8 leaves this open; the scan answers it. Cross-checked here
# against next-generation matrices built directly, sharing no code with
# ngm_switched (which routes through the f_T/f_D/f_C factorisation).
from theory import next_gen_matrix, spectral_radius, degree_moments
mean_k, cross = degree_moments(P, 2)

def rho_mf_direct(m, lam):
    X = cross / mean_k[:, None]
    return spectral_radius(np.array([[(m - 1) * lam * X[a, b] for b in range(2)]
                                     for a in range(2)]))

for m, lam, want_flip in [(8, 0.2, True), (12, 0.2, True), (16, 0.15, True),
                          (4, 0.2, False), (8, 0.02, False), (8, 1.0, False)]:
    r_ex = spectral_radius(next_gen_matrix(P, (m, m), lam))
    r_mf = rho_mf_direct(m, lam)
    allpass &= ok(f"sign of net deviation at m={m}, lam={lam} "
                  f"({'MF under' if want_flip else 'MF over'})",
                  (r_ex > r_mf) == want_flip,
                  f"rho={r_ex:.4f} vs rho_MF={r_mf:.4f}")
    # the factor path must agree with the direct matrices
    allpass &= ok("   factor path == direct matrices",
                  abs(r_ex / r_mf - factors(P, (m, m), lam)[3]) < 1e-10)

first = None
for m in range(2, 25):
    ls = np.geomspace(1e-3, 20, 2000)
    if any(factors(P, (m, m), l)[3] > 1 for l in ls):
        first = m
        break
allpass &= ok("first group size with a sign flip is m=8", first == 8, f"m={first}")

# at the threshold itself the flip never occurs in the explored range
no_flip_at_lc = all(factors(P, (m, m), lambda_c_rho(P, (m, m)))[3] < 1
                    for m in range(2, 21))
allpass &= ok("no sign flip AT the threshold for m<=20 (so lam_c^MF < lam_c)",
              no_flip_at_lc)

print("\n=== ALL PASS ===" if allpass else "\n!!! SOME FAILURES !!!")
