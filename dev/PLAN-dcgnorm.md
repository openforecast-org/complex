# Plan: complex generalised normal distribution (`dcgnorm`) in `complex`

Branch: `dcgnorm` (from `master` at `a8b660c`, v1.0.2). Prototype with numerical checks:
`dev/prototype-dcgnorm.R`. `dev/` is in `.Rbuildignore`.

Motivation: residuals of wind data have much heavier tails than normal (see the cARIMA wind
project). The complex normal likelihood in `clm()` is the only likelihood available. Goal: add the
complex generalised normal (CGN) with shape β estimated, and only ONE extra estimated parameter,
by concentrating the scale and pseudo-scale out analytically.

## 1. The distribution

Parameters: location μ (complex), scale α² > 0 (real), pseudo-scale ψ² (complex, |ψ²| < α²),
shape β > 0. For e = z − μ:

    Q(e)  = 2 (α²|e|² − Re(Conj(ψ²) e²)) / (α⁴ − |ψ²|²)
    f(z)  = 2β / (π Γ(1/β) 2^(1/β) sqrt(α⁴ − |ψ²|²)) · exp(−½ Q^β)

- β = 1: complex normal with σ² = α², ς² = ψ² (checked in the prototype: max log-density difference from `dcnorm` is 2e-16).
- β < 1 heavier tails (β = 0.5 Laplace-like), β → ∞ uniform on an ellipse.
- Moments: variance σ² = k(β) α², pseudo-variance ς² = k(β) ψ², where
  k(β) = 2^(1/β) Γ(2/β) / (2 Γ(1/β)). Hence ς²/σ² = ψ²/α² for all β.
- Radial law: Q^β ~ Gamma(shape = 1/β, scale = 2). Used for simulation and diagnostics.
- Equivalent to the bivariate exponential power distribution (Gómez et al., 1998) and to the
  complex generalized Gaussian of Novey, Adalı & Roy (2010, IEEE TSP 58(3)) under the
  `cnorm.R` mapping Σ = ½[[α²+Re ψ², Im ψ²],[Im ψ², α²−Re ψ²]].

**Parameterisation decision (to confirm):** user-facing functions take `sigma2`, `varsigma2`
(variance and pseudo-variance, as in `dcnorm`) plus `shape`, and convert internally with k(β).
This keeps `sigma.clm()` meaningful and `dcgnorm(shape=1)` == `dcnorm`. Alternative: expose
`scale2`/`pseudoscale2` directly. Recommendation: variance parameterisation.

## 2. Concentrated likelihood (closed form)

Write ψ² = ρ α², ρ the complex circularity coefficient (|ρ| < 1).

1. ρ̂ = mean(e²) / mean(|e|²). Moment estimator, consistent for every β, exact MLE when β = 1.
2. Given ρ̂, qₜ = 2(|eₜ|² − Re(Conj(ρ̂) eₜ²)) / (1 − |ρ̂|²), and the exact MLE of the scale is
   α̂² = [ (β / (2T)) Σ qₜ^β ]^(1/β),  ψ̂² = ρ̂ α̂².
3. Plugging in gives, since Σ (qₜ/α̂²)^β = 2T/β,
   ℓ(B, β) = T [ log(2β) − log π − lgamma(1/β) − log(2)/β − log α̂² − ½ log(1 − |ρ̂|²) − 1/β ].
   With β = 1 this equals the current `clm` cost, T(log 2π + 1 + ½ log det Σ̂).

Prototype check (T = 20,000, β = 0.5, α² = 2, ψ² = 0.6+0.8i): closed form ℓ = −101981.6 vs
full numeric MLE −101981.5; the β that maximises the closed-form likelihood is 0.50.

Optional exact-MLE refinement (argument `refine`, default 0): n fixed-point steps
wₜ = Qₜ^(β−1), α² = (β/T) Σ wₜ|eₜ|², ψ² = (β/T) Σ wₜ eₜ², then recompute ℓ with the full
density. Still closed-form arithmetic; the optimiser sees only β. Convergence is proven for
β ≤ 1 (Pascal et al., 2013); for β > 1 cap the iterations and fall back to step 2.

## 3. Implementation phases

Each phase ends with `R CMD build` + `R CMD check --no-manual` clean and its tests passing.
Commit per phase.

### Phase 1 — distribution functions (`R/cgnorm.R`, new)
- `dcgnorm(q, mu=0, sigma2=1, varsigma2=0, shape=1, log=FALSE)`. Default shape 1 = normal.
  In greybox `dgnorm` shape 2 = normal (exponent on |x|, not on the quadratic form); document
  that shape here equals greybox shape / 2.
- `rcgnorm(n, mu, sigma2, varsigma2, shape)` via the Gamma radial law and a uniform angle.
- Internal helper `cgnormScale(sigma2, varsigma2, shape)` returning (α², ψ²), and its inverse.
- Validate: sigma2 > 0, |varsigma2| < sigma2, shape > 0.
- No `p`/`q` functions in this release (no closed form; `qcnorm` semantics are unusual anyway).
- Roxygen docs in the `cnormal` style; export; NAMESPACE via roxygen.

### Phase 2 — `clm()` estimation (`R/clm.R`)
- New argument `distribution = c("dcnorm", "dcgnorm")`, used only when `loss = "likelihood"`.
  Default keeps current behaviour exactly.
- Shape: pass via `...` as `shape` (NULL = estimate, number = fix), mirroring greybox `alm()`
  handling of `other`. Returned in `object$other$shape`.
- `fitter()` (≈ lines 195–247): for `dcgnorm`, compute ρ̂, α̂², ψ̂² as in section 2 and return
  `scale` as the 2×2 covariance (via k(β)), so downstream code keeps its meaning.
- Parameter vector: append ONE real parameter log β at the end of `BReal` in `estimator()`
  (≈ lines 307–341) and strip it BEFORE the real/imaginary split in `fitter()` (lines 197–201),
  which currently assumes every entry is half of a complex coefficient. Update `maxeval`.
- `CF()` (≈ lines 259–305): `dcgnorm` branch = −ℓ(B, β) from section 2.
- Starting values: B from the current `dcnorm` fit logic, log β = 0.
- `nParam` (line 925): add 0.5 when β is estimated (package counts in complex units).
- `logLik` stored as −CF; `AIC`/`AICc` then work via `nparam`/`logLik` methods.

### Phase 3 — inference (`vcov.clm`, ≈ lines 1041–1215)
- Hessian (line 901) currently over B only. For `dcgnorm` with estimated β compute over (B, log β)
  and take the B block of the inverse. The refit call in `vcov.clm` must pass `distribution`
  and `shape`.
- Delta method for the SE of β if reported in `summary`.

### Phase 4 — prediction (`predict.clm`, ≈ lines 1368–1600)
- Point forecasts unchanged.
- Intervals for `dcgnorm`: simulate. Draw ε from `rcgnorm` with fitted (σ², ς², β), push through
  the ARIMA recursion for h steps, take quantiles of Re and Im separately (matching the current
  per-part intervals). `nsim` argument, default 10000. Return the simulated paths optionally
  (needed for energy scores in the wind paper).
- Leave the `dcnorm` interval code as is (only warn as now).

### Phase 5 — output and diagnostics
- `print.summary.clm`: print distribution and shape (with SE when estimated).
- `plot.clm`: the greybox normal Q–Q plot is wrong for `dcgnorm`; add a Q–Q plot of Qₜ^β against
  Gamma(1/β, 2) (a new `which` value), or at least skip the normal Q–Q plot.
- NEWS entry, DESCRIPTION version bump (1.1.0), README mention.

## 4. Tests (`tests/testthat/`, none exist yet; `testthat` already in Suggests)

1. `dcgnorm(shape=1)` == `dcnorm` (tolerance 1e-12) for several (sigma2, varsigma2).
2. Density integrates to 1 on a grid for shape ∈ {0.3, 0.5, 1, 2, 5}.
3. `rcgnorm`: sample σ², ς² within 3 SE of targets; KS test of Q^β against Gamma(1/β, 2) p > 0.001
   (fixed seed).
4. Closed-form ℓ equals −Σ log `dcgnorm` evaluated at (α̂², ψ̂²) (identity check of section 2).
5. `clm(y~x, distribution="dcgnorm", shape=1)` reproduces `clm(y~x)` coefficients (1e-6) and logLik.
6. cAR(1) simulated with β = 0.5, T = 2000: β̂ in [0.4, 0.6], AICc(dcgnorm) < AICc(dcnorm).
7. Regression: every existing example in the docs gives unchanged results with the default
   `distribution`.

## 5. Open questions for the author

1. Variance parameterisation (recommended) or scale parameterisation for `dcgnorm`?
2. `distribution` as a new argument (recommended) or extend `loss` values (e.g. `"likelihood-cgn"`)?
3. Parameter counting: `nparam` counts complex units (one complex coefficient = 1, Σ = 1.5).
   AIC/AICc of `clm` therefore penalise half the real-parameter count. Keep the convention (and
   document it) or switch to real counts? This affects comparing cARIMA with real VAR by IC.
4. Ship `refine` (exact MLE of ρ) in the first version or later?
5. Name: `dcgnorm` (parallel to greybox `dgnorm`) — confirm.
