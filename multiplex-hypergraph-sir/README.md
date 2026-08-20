# Multiplex-hypergraph SIR — group closure & the subcritical determination of λc

Reference implementation of the group-closure theory for SIR on multiplex
hypergraphs, and the determination of the epidemic threshold **λc by
extrapolation from subcritical final states** (research outline §3.1, Figure 2).

The threshold is measured a *second, independent way* from the eigenvalue
condition ρ(N)=1: the nonlinear closure is integrated to its subcritical final
state, and the reciprocal outbreak size ε/R(∞) is extrapolated linearly to its
zero-crossing. The two paths share no code and agree to **6.5×10⁻⁴ (≤ 1.51σ)**
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
| `simulate.py`       | exact Gillespie SIR on a configuration-model multiplex hypergraph (real loops, finite N) — the microscopic ground truth |
| `figure2_evolution.py` + `figure2_evolution.{pdf,png}` | time-evolution figure: group closure vs exact Gillespie simulation for S(t), I(t) (Fenwick-tree Gillespie; N=6000, 400 runs) |
| `figure_lambda_c.py`| draws the λc figure from `figure_data.json` |
| `figure3_lambda_c.pdf/.png` | the λc figure (extrapolation vs ρ(N)=1) |
| `figure4_rho12.py` + `figure4_rho12.{pdf,png}` | Figure 4: λc vs inter-layer participation correlation ρ12 — exact ρ(N)=1 curve vs bootstrap box plots of the Gillespie threshold; sim data cached in `figure4_data.json` |

## Reproduce

```bash
pip install -r requirements.txt
python3 verify.py         # analytic anchors           -> ALL PASS
python3 verify_ode.py     # closure self-checks        -> ALL PASS
python3 datacheck.py      # independent data checks     -> agreement < 1e-3
python3 recheck_mc.py     # single-group Gillespie MC    -> C matches within noise
python3 figure2_evolution.py  # closure vs Gillespie, S(t)/I(t) (slow)
python3 lambda_c_extrap.py# recompute figure_data.json
python3 figure_lambda_c.py# redraw the λc figure
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
with **max relative deviation 6.5×10⁻⁴**, every configuration consistent with
ρ(N)=1 within **1.51σ** of the regression uncertainty (the eight pulls are
+1.51, +1.38, +1.26, +1.03, +0.76, +0.57, +0.55, +0.47). The residual is the known
finite-window curvature of the linear law (1/χ is strictly linear only near λc);
its sign flips with window placement and its size is comparable to σ, so it
cannot be distinguished from a true deviation.

## Time-evolution validation

Beyond the threshold, the group closure reproduces the full transient. At
M=2, m=(3,4), N=6000, λ=1.6λc, ε=0.02, the closure S(t)/I(t) (Eq. 2.5 plus
one quadrature Ṙ=μI to split I from R) sits on top of an exact continuous-time
Gillespie simulation (`figure2_evolution.py`, Fenwick-tree group sampling,
400 independent runs): max deviation 0.0019 in S and 0.0006 in I — both under a
line width — with matching peak height (~0.3%) and peak time (t≈4.4). The
residual is a finite-size effect that vanishes as N grows.
