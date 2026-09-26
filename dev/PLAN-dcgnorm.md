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

**Parameterisation (DECIDED):** scale parameterisation. Functions take `scale` (α², real) and
`pseudoscale` (ψ², complex) plus `shape` (β). The closed form is simplest in the scale (section 2);
the variance form only adds the factor k(β). At shape = 1, scale = σ² and pseudoscale = ς², so
`dcgnorm(q, mu, scale=s, pseudoscale=p, shape=1)` == `dcnorm(q, mu, sigma2=s, varsigma2=p)`.
`sigma.clm()` keeps returning the covariance of residuals; `object$scale` holds (α², ψ²) for dcgnorm.

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

Exact MLE of ρ (fixed point wₜ = Qₜ^(β−1), α² = (β/T) Σ wₜ|eₜ|², ψ² = (β/T) Σ wₜ eₜ²) is NOT
implemented (DECIDED). Document in `?clm` and `?dcgnorm` that the scale/pseudo-scale estimates are
approximate ML for β ≠ 1: ρ is the moment estimator, the scale is exact ML given ρ.

## Status (2026-09-25)

Phase 0 done (IC fix). Phases 1 and 2 done: `R/cgnorm.R` (dcgnorm, rcgnorm, internal
cgnormConcentrated), clm(distribution=), shape estimated as log-shape appended to the optimiser
vector, nparam +0.5, AICc/BICc m excludes the shape, vcov refit fixes the shape (Hessian for B
given shape, Phase 3 still to do properly), summary prints the shape. Tests in
tests/testthat/test-cgnorm.R. Remaining: Phase 3 (joint Hessian), Phase 4 (simulated intervals),
Phase 5 (diagnostic plot), and clm's poor default starting values (found in the wind experiment:
nested bigger models can end with lower logLik).

## 3. Implementation phases

Each phase ends with `R CMD build` + `R CMD check --no-manual` clean and its tests passing.
Commit per phase.

### Phase 1 — distribution functions (`R/cgnorm.R`, new)
- `dcgnorm(q, mu=0, scale=1, pseudoscale=0, shape=1, log=FALSE)`. Default shape 1 = normal.
  In greybox `dgnorm` shape 2 = normal (exponent on |x|, not on the quadratic form); document
  that shape here equals greybox shape / 2.
- `rcgnorm(n, mu, scale, pseudoscale, shape)` via the Gamma radial law and a uniform angle.
- Exported helper or documented formula k(β) to convert scale to covariance.
- Validate: scale > 0, |pseudoscale| < scale, shape > 0.
- No `p`/`q` functions in this release (no closed form; `qcnorm` semantics are unusual anyway).
- Roxygen docs in the `cnormal` style; export; NAMESPACE via roxygen.

### Phase 2 — `clm()` estimation (`R/clm.R`)
- New argument `distribution = c("dcnorm", "dcgnorm")` (DECIDED), as in greybox `alm()`: the
  distribution defines the likelihood used when `loss = "likelihood"`. Default keeps current
  behaviour exactly.
- Shape: as `alm()` does for `dgnorm`: user may fix it via `...` (`shape=`), otherwise estimated.
  Returned in `object$other$shape`.
- `fitter()` (≈ lines 195–247): for `dcgnorm`, compute ρ̂, α̂², ψ̂² as in section 2 and return
  them as `scale`. Check every downstream use of `object$scale` (vcov refit passes it back in).
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

## 5. Information criteria: parameter count (found while planning, affects master too)

Check (cAR(1) example from the monograph, T = 80):
- `logLik.clm` is the FULL bivariate log-likelihood, −T(log 2π + 1 + ½ log det Σ̂) (−389.60).
- `nparam.clm` returns `object$rank` = `nVariables + 3/2` = 3.5, i.e. real parameters / 2
  (2 intercept + 2 AR + 3 Σ = 7 real). Commented-out code in `nparam.clm` says
  "Divide by two, because we calculate df per series".
- `AIC`/`BIC` (stats, via `logLik` df attr) and `AICc`/`BICc` (greybox `.default`, via `nparam`)
  all use 3.5. AIC = −2ℓ + 7 = 786.21; with all 7 real parameters it is −2ℓ + 14 = 793.21.
- For comparison greybox's own `AICc.varest` counts ALL real parameters of a VAR in the penalty,
  K·m + K(K+1)/2 (coefficients of all equations plus Σ), and uses the per-equation count m only in
  the small-sample denominator (Bedrick & Tsai 1994): −2ℓ + 2T(Km + K(K+1)/2)/(T − m − K − 1).
  (`vars::logLik.varest` reports df = 6 for this VAR(1), omitting Σ; do not compare to it blindly.)

Conclusion: the per-series count (3.5) is right for `df.residual` (each series has T obs and
3.5 parameters per series), but with a joint bivariate log-likelihood the AIC/BIC penalty needs
the total real count (7). The current values under-penalise by half, so AIC-based order
selection favours bigger cARIMA models, and IC comparison with real-valued VAR favours cARIMA.

Proposed fix (Phase 0, separate commit, NEWS entry as a behaviour change):
- Keep `nparam.clm` / `df.residual` per series (no change to t-based inference).
- `logLik.clm`: `df = 2*nparam(object)` so stats `AIC`/`BIC` use the total real count.
- New `AICc.clm`, `BICc.clm` in the `AICc.varest` style with K = 2, k = 2*nparam(object) and
  m = nparam(object) − (Σ and shape share) in the denominator. The exact small-sample
  correction for the complex (restricted) model is not derived anywhere; decide the form of m.
- dcgnorm: β adds 1 real parameter (0.5 in the per-series count).
- Re-run the monograph examples: selected orders may change (BJsales chose cARIMA(3,1,1)).

## 6. Decisions (2026-09-25)

1. Scale parameterisation (`scale`, `pseudoscale`).
2. `distribution` argument, as in `alm()`.
3. IC count: approved and implemented as proposed in section 5 (Phase 0 done, per-series nparam kept,
   legion-style AICc/BICc). BJsales check: top-3 orders by AICc and BIC unchanged before/after
   (cARIMA(3,1,2) first. The monograph reported (3,1,1); the old formula on the same fits also
   gives (3,1,2), so the difference predates this fix. Cause not investigated.)
4. No exact-MLE refinement; document the approximation.
5. Name `dcgnorm` / `rcgnorm`.
