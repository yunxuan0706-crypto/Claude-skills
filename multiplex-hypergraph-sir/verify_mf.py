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
# includes ASYMMETRIC layers: the per-layer definition of f_T/f_C used to pass
# this test for m1==m2 while giving the wrong sign at e.g. m=(2,10)
for m in [(3, 3), (5, 5), (2, 10), (3, 12), (2, 8), (4, 9)]:
    for lam in (0.05, 0.2, 1.0):
        fT, fD, fC, net = factors(P, m, lam)
        prod = fT * fD * fC
        allpass &= ok(f"factors multiply to net (m={m}, lam={lam})",
                      abs(prod - net) < 1e-10,
                      f"f_T={fT:.4f} f_D={fD:.4f} f_C={fC:.4f} -> {prod:.6f} vs {net:.6f}")
        allpass &= ok(f"  directions  f_T<1, f_D<1, f_C>1 (m={m}, lam={lam})",
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

def first_flip(Pd):
    ls = np.geomspace(1e-3, 20, 1200)
    for mm in range(2, 41):
        if any(factors(Pd, (mm, mm), l)[3] > 1 for l in ls):
            return mm
    return None

allpass &= ok("first group size with a sign flip is m=8 for THIS config",
              first_flip(P) == 8, f"m={first_flip(P)}")

# ---- 8. how far do those statements generalise? ---------------------------
# The onset of the flip is config dependent (it is set by f_D, so more
# heterogeneous degrees flip earlier); the absence of a flip AT the threshold
# is not. Both are reported honestly rather than assumed.
OTHER = {
    "k in {3,5} independent": {(3, 3): .25, (3, 5): .25, (5, 3): .25, (5, 5): .25},
    "3-regular":              {(3, 3): 1.0},
    "10-regular":             {(10, 10): 1.0},
    "k in {2,20} heavy tail": {(2, 2): .9, (20, 20): .1},
    "k in {1,3} sparse":      {(1, 1): .5, (3, 3): .5},
}
onsets = {"k in {3,5} corr": first_flip(P)}
for nm, Pd in OTHER.items():
    onsets[nm] = first_flip(Pd)
allpass &= ok("a sign flip exists for every distribution tested",
              all(v is not None for v in onsets.values()),
              "onsets " + ", ".join(f"{k}:{v}" for k, v in onsets.items()))
allpass &= ok("flip onset lies in m = 6..9 across those distributions",
              all(6 <= v <= 9 for v in onsets.values()))

# f_D is NOT (X-1)/X. Every X here has positive entries, so subtracting the
# identity shifts all eigenvalues by -1 and f_D = 1 - 1/rho(X) exactly; for
# M=2 uniform-X configs that is (2X-1)/(2X), not (X-1)/X. The wrong form is
# checked too, so a regression to it cannot pass silently.
print("\nf_D vs the degree ratio, and the onset it sets:")
pairs = []
for nm, Pd in list(OTHER.items()) + [("k in {3,5} corr", P)]:
    mk, cr = degree_moments(Pd, 2)
    Xm = cr / mk[:, None]
    rX = spectral_radius(Xm)
    fD = factors(Pd, (5, 5), 0.2)[1]
    allpass &= ok(f"f_D = 1 - 1/rho(X)  [{nm}]", abs(fD - (1 - 1 / rX)) < 1e-12,
                  f"f_D={fD:.6f}  1-1/rho(X)={1 - 1/rX:.6f}  "
                  f"(X-1)/X={(Xm[0,0]-1)/Xm[0,0]:.6f}")
    allpass &= ok(f"   and f_D != (X-1)/X  [{nm}]",
                  abs(fD - (Xm[0, 0] - 1) / Xm[0, 0]) > 1e-3)
    pairs.append((rX, onsets[nm], nm))
pairs.sort()
allpass &= ok("onset is monotone non-increasing in rho(X)",
              all(pairs[i][1] >= pairs[i + 1][1] for i in range(len(pairs) - 1)),
              "  ".join(f"{nm}:rho={r:.2f}->m={o}" for r, o, nm in pairs))

for nm, Pd in list(OTHER.items()) + [("k in {3,5} corr", P)]:
    no_flip = all(factors(Pd, (m, m), lambda_c_rho(Pd, (m, m)))[3] < 1
                  for m in range(2, 31))
    allpass &= ok(f"no flip AT threshold, m<=30  [{nm}]", no_flip)

# ---- 9. analytic identity at m = 2 ----------------------------------------
# For a config whose X_ab = <k^a k^b>/<k^a> is the same for every (a,b), the
# m=2 threshold ratio collapses to f_D exactly:
#   rho_ex = T(2X-1), rho_MF = 2X*lam  =>  lam_c = 1/(2X-2), lam_c^MF = 1/(2X)
#   =>  lam_c^MF / lam_c = (X-1)/X.
# (It does NOT hold when the layers are independent, since then X_aa != X_ab.)
for nm, Pd in [("k in {3,5} corr", P), ("3-regular", {(3, 3): 1.0}),
               ("10-regular", {(10, 10): 1.0}),
               ("k in {2,20}", {(2, 2): .9, (20, 20): .1})]:
    mk, cr = degree_moments(Pd, 2)
    X = cr[0, 0] / mk[0]
    r = lambda_c_mf(Pd, (2, 2)) / lambda_c_rho(Pd, (2, 2))
    allpass &= ok(f"m=2 identity lam_c^MF/lam_c = (X-1)/X  [{nm}]",
                  abs(r - (X - 1) / X) < 1e-9,
                  f"{r:.10f} vs {(X-1)/X:.10f}")

# ---- 10. mean-field ODE closes the loop around its own threshold -----------
# The Jacobian check above is a linearisation; integrate the mean-field ODE
# itself and confirm it is sub/supercritical on either side of lam_c^MF.
from scipy.integrate import solve_ivp
for m in (3, 6):
    lc = lambda_c_mf(P, (m, m))
    for frac, grows in ((0.85, False), (1.15, True)):
        mf = MeanFieldODE(P, (m, m), frac * lc)
        K = mf.K
        y0 = np.concatenate([np.full(K, 1 - 1e-6), np.full(K, 1e-6)])
        sol = solve_ivp(mf.rhs, [0, 400], y0, method="LSODA",
                        rtol=1e-10, atol=1e-14)
        final_i = sol.y[K:, -1].max()
        peak_i = sol.y[K:, :].max()
        took_off = peak_i > 1e-3
        allpass &= ok(f"MF ODE {'super' if grows else 'sub'}critical at "
                      f"{frac}*lam_c^MF (m={m})", took_off == grows,
                      f"peak i={peak_i:.2e}, final i={final_i:.2e}")

print("\n=== ALL PASS ===" if allpass else "\n!!! SOME FAILURES !!!")
