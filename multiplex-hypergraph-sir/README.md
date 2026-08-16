# Multiplex-hypergraph SIR — group closure & the subcritical determination of λc

Reference implementation of the group-closure theory for SIR on multiplex
hypergraphs, and the determination of the epidemic threshold **λc by
extrapolation from subcritical final states** (research outline §3.1, Figure 2).

The threshold is measured a *second, independent way* from the eigenvalue
condition ρ(N)=1: the nonlinear closure is integrated to its subcritical final
state, and the reciprocal outbreak size ε/R(∞) is extrapolated linearly to its
zero-crossing. The two paths share no code and agree to **6×10⁻⁴ (≤ 1.5σ)**
across eight configurations.

## Files

| file | what it is |
|------|------------|
| `theory.py`         | core theory: cascade `C(m,λ,θ)` (2.8–2.9), next-generation matrix `N` (2.10), `ρ(N)=1` threshold, and the group-closure ODE (2.2, 2.4–2.6) |
| `lambda_c_extrap.py`| the deliverable: subcritical-final-state extrapolation of λc, with ε→0 Richardson; run as a script to regenerate `figure_data.json` |
| `verify.py`         | reproduces the outline's analytic anchors (cascade values, Anchor A ρ=1.013, degenerate λc=1/3, θ≥2 drop-out) |
| `verify_ode.py`     | closure self-checks: sum rule (2.6), identity (2.7), DFE Jacobian abscissa = 0 at λc |
| `datacheck.py`      | independent cross-checks of the final-state data (branching formula, ε→0, t_max, sanity) |
| `recheck_mc.py`     | direct single-group Gillespie MC for the cascade `C`; ODE behaviour across the transition; `ρ(λ)` monotonicity |
| `meanfield.py`      | degree mean-field closure (2.12–2.13) and the switchable next-generation matrix that turns each of the three deviations on/off |
| `verify_mf.py`      | mean-field checks: the λc^MF/λc anchor, MF-matrix vs MF-ODE Jacobian, factor directions, and the sign-flip scan |
| `simulate.py`       | exact Gillespie SIR on a configuration-model multiplex hypergraph (real loops, finite N) — the ground truth for Figure 3 |
| `verify_sim.py`     | simulation validation of Figure 3: χ_sim vs exact N vs mean field, cascade C sampled at the flip-region parameters, and the supercriticality of the flip |
| `figure3_meanfield.py` + `figure3_meanfield.{pdf,png}` | Figure 3: the three factors, the net deviation on the (m,λ) plane, and the threshold consequence |
| `figure_lambda_c.py`| draws Figure 2 from `figure_data.json` |
| `figure2_lambda_c.pdf/.png` | Figure 2 |

## Reproduce

```bash
pip install -r requirements.txt
python3 verify.py         # analytic anchors           -> ALL PASS
python3 verify_ode.py     # closure self-checks        -> ALL PASS
python3 datacheck.py      # independent data checks     -> agreement < 1e-3
python3 recheck_mc.py     # single-group Gillespie MC    -> C matches within noise
python3 verify_mf.py      # mean-field checks            -> ALL PASS
python3 verify_sim.py     # Gillespie ground truth       -> ALL PASS (slow)
python3 figure3_meanfield.py  # redraw Figure 3
python3 lambda_c_extrap.py# recompute figure_data.json
python3 figure_lambda_c.py# redraw the figure
```

## Method (Figure 2)

In the subcritical regime the mean outbreak size from a small seed ε is
`χ(λ) = 1ᵀ(I−N)⁻¹v₀`, which diverges as `(1−ρ)⁻¹` as `ρ(N)→1⁻`. Hence its
reciprocal `ε/R(∞) = 1/χ` falls **linearly** to zero and its zero-crossing is
λc. `R(∞) = 1 − (1−ε)Ψ(Φ(∞))` is the final size of the nonlinear closure
(Eq. 2.5) integrated to steady state — an independent numerical path from the
ρ(N)=1 eigenvalue computation. The finite-seed bias (a small seed still yields a
nonzero outbreak fraction exactly at λc, pushing the zero-crossing high) is
removed by Richardson extrapolation ε→0, since R(∞)/ε is linear in ε.

## Verification summary

**Analytic anchors** (`verify.py`): C(2,λ,1)=λ/(1+λ) exactly; C(m,1,1) exceeds
(m−1)T by 11/22/31% for m=3,4,5; C=0 for θ≥2; **Anchor A** P={(2,2),(3,3)},
m=(3,3), λ=0.13 → N₁₁=0.386, ρ=1.013; **degenerate** 5-regular pairwise →
λc=1/3 to 10 digits; **Anchor B** a θ=2 channel drops out of ρ(N)=1 exactly
(the outline's 0.410619 is the finite-t_max ODE-bisection value, biased +2.51%).

**Closure self-checks** (`verify_ode.py`): sum rule (2.6) residual ~1e-17;
identity (2.7) residual ~1e-12; the disease-free Jacobian's leading eigenvalue
crosses zero exactly at λc (machine precision) — the closure and N share a
threshold.

**Independent data checks** (`datacheck.py`): the closure final states match the
branching formula `1 + 1ᵀ(I−Nᵀ)⁻¹g`, `g_b=⟨k_b⟩C_b`, converging as ε→0
(3.7e-3 → 1.1e-5 for ε=2e-4 → 1e-5); final states are t_max-converged
(χ identical over t_max∈[2e4,1.2e5], activated-mass residual ≤ 1e-20).

**Monte-Carlo re-check** (`recheck_mc.py`): a direct single-group Gillespie
simulation reproduces the cascade `C(m,λ,θ)` within statistical noise (≤ ~1σ;
the script samples 2×10⁵ per point, and a 10⁶ × 4-seed re-check of the m=2 case
gave 0.64–1.12σ), incl. `C=λ/(1+λ)` for m=2 and `C=0` for θ≥2); the closure `R(∞)`
stays O(ε) below λc and lifts to an O(1) epidemic above it; `ρ(λ)` is strictly
monotone (bisection well-posed).

## Result

Across eight configurations spanning λc ∈ [0.08, 0.40] (single-layer higher-order,
pairwise-degenerate, and two-layer), the subcritical extrapolation recovers λc
with **max relative deviation 6×10⁻⁴**, every configuration consistent with
ρ(N)=1 within **1.5σ** of the regression uncertainty. The residual is the known
finite-window curvature of the linear law (1/χ is strictly linear only near λc);
its sign flips with window placement and its size is comparable to σ, so it
cannot be distinguished from a true deviation.

## Mean-field deviations, term by term (Figure 3)

The degree mean-field next-generation matrix `N^MF_ab = (m_b−1)λ_b⟨k^a k^b⟩/⟨k^a⟩`
differs from the exact `N_ab` in exactly three places, with two competing
directions. Written as multiplicative factors on the spectral radius:

| switch | what mean field gets wrong | factor | direction |
|---|---|---|---|
| **T** | uses λ instead of `T=λ/(1+λ)` — ignores that a group is used up on a member once it transmits | `f_T = 1/(1+λ) < 1` | MF **over**estimates ρ |
| **D** | omits the excess subtraction `−δ_ab` — lets infection turn back along the group it arrived by | `f_D < 1` (0.882 here) | MF **over**estimates ρ |
| **C** | uses `(m−1)T` instead of the cascade `C` — cannot see that intra-group infections lengthen the active period | `f_C = C/[(m−1)T] > 1` | MF **under**estimates ρ |

`ngm_switched(..., s_T, s_D, s_C)` turns each on independently; `(1,1,1)` is
exact and `(0,0,0)` is mean field, verified to 1e-9.

**Result.** With `P={(3,3),(5,5)}` (k∈{3,5}, two layers, m₁=m₂=m):

- λc^MF/λc = 0.765, 0.849, 0.879, 0.903, 0.913 for m = 2,3,4,6,8, rising
  monotonically to 0.93 at m=16 — mean field **always underestimates λc**, and
  the gap narrows with m without changing sign. At m=3 the term-by-term
  thresholds are ×1.063 (T only), ×1.133 (D only), ×0.976 (C only) relative to
  mean field, so D is the largest single term and C opposes the other two.
- **The net deviation of ρ does change sign** on the (m,λ) plane. Above a
  critical group size, at moderate λ, the cascade term wins and mean field
  *under*estimates ρ. For the figure's configuration the onset is m = 8, with
  window λ∈[0.086, 0.394], widening with m (max `ρ/ρ^MF` = 1.58 at m=16).
  The outline lists this as an open question; the scan answers it in the
  affirmative, and `verify_mf.py` confirms the sign from next-generation
  matrices built directly, sharing no code with the factorisation.
- **The onset is configuration dependent, the existence of the flip is not.**
  It is set by `f_D = (X−1)/X` with `X = ⟨k^a k^b⟩/⟨k^a⟩`, so more heterogeneous
  degrees flip earlier: across six distributions the onset spans m = 6–9
  (m=6 for 10-regular and for a heavy tail k∈{2,20}; m=9 for 3-regular and
  k∈{1,3}), and a flip exists in every one of them.
- **But the flip never reaches the threshold.** λc lies far below the flip
  window for every m tested (λc = 0.018 vs window [0.086, 0.394] at m=8), so
  λc^MF < λc always — verified for all six distributions up to m = 30. The
  window's lower edge, in units of λc, falls with m (4.7 → 2.1 between m=8 and
  m=50) but stays clear of 1. The
  threshold statement and the reproduction-number statement are genuinely
  different, and only the latter changes sign.
- At m = 2 the threshold ratio collapses to an exact identity,
  `λc^MF/λc = (X−1)/X`, for any configuration with a uniform `X_ab`
  (0.667 for 3-regular, 0.765 here, 0.900 for 10-regular; verified to 1e-9).
  It does not hold when the layers are independent, since then `X_aa ≠ X_ab`.

The comparison is for θ = 1 throughout: the mean-field closure (2.12) takes
`Θ_a = 1−(1−φ_a)^{m_a−1}`, which presumes a single infected member activates a
group. For θ ≥ 2 the exact side has `C = 0` and the two closures are not
comparable in this form.

### Simulation ground truth for Figure 3

Figure 3 compares two closures, so on its own it is theory against theory. The
missing ground truth is supplied in `verify_sim.py`, in two parts dictated by a
structural fact: **the sign-flip region is always supercritical** — its lower
edge stays above λc at every m tested — the ratio decreases with m (about 4.7 at
m=8, 2.6 at m=16, 2.1 by m=50) but never approaches 1 — so it cannot be
reached by a subcritical outbreak-size measurement.

1. **The full next-generation bookkeeping is validated where it can be.** Exact
   Gillespie SIR on a configuration-model multiplex hypergraph (N = 40 000,
   3 000 seeds spread over 3 independent graphs, real loops and finite-N
   effects, no closure) at
   λ = 0.6 λc and 0.8 λc for m = 3, 5, 8, 12. The simulated mean outbreak size
   agrees with `χ = 1ᵀ(I−N)⁻¹v₀` at **|z| ≤ 1.3σ on all eight points**, while
   the mean-field prediction is excluded by **3.5σ to 42.7σ**. This tests T,
   the excess subtraction and the cascade together. (At N = 12 000 one point
   sat at −3.1σ; raising N moved it to −0.4σ, confirming finite-N, not a
   discrepancy.)
2. **The term that drives the flip is validated inside the flip region.** The
   cascade `C(m,λ)` is sampled by direct single-group Gillespie at the
   flip-region parameters themselves — (m,λ) = (8, 0.15), (8, 0.30), (10, 0.20),
   (12, 0.30), (16, 0.15), all with net > 1 — and matches the recursion to
   ≤ 1.1σ at 2×10⁵ realisations.

So every ingredient the flip is built from is independently validated; the flip
itself is a prediction of that validated matrix in a regime where the epidemic
has already taken off, and is labelled as such.