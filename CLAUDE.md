# complex (R package)

Complex-valued econometrics and forecasting: complex linear models, complex ARIMA (cARIMA), complex
distributions, complex correlations. Companion to Svetunkov & Svetunkov (2024), *Complex-Valued
Econometrics with Examples in R* (Springer); Chapter 6 covers the dynamic models.
Author/maintainer: Ivan Svetunkov. Depends on greybox and legion (which loads smooth).

## Repository

- `git@github.com:openforecast-org/complex.git` (moved from config-i1/complex). Always use ssh.
  There is no ssh in the Claude container, so the user does push/fetch.
- Branch `dcgnorm`: development of v1.0.3 (not yet merged to `master`). Commit as
  `Ivan Svetunkov <ivan@svetunkov.com>` (`git -c user.name=... -c user.email=...`).
- Used by the wind forecasting paper project: `~/R/Projects/Experiments/cARIMA`
  (git@github.com:config-i1/complex-wind.git), which has its own CLAUDE.md.

## Layout

| file | content |
|---|---|
| `R/clm.R` | `clm()` (complex linear model and cARIMA(X)) and its methods: `logLik`, `nparam`, `AICc/BICc`, `clmIC()`, `vcov`, `confint`, `summary`, `predict`, `plot` |
| `R/carima.R` | `carima()` pure cARIMA wrapper of `clm()` (msarima-like interface), `forecast.carima`, `carimaFinalise()` |
| `R/autoCarima.R` | `auto.carima()` order selection (stepwise / full / fast) |
| `R/hannanRissanen.R` | Hannan-Rissanen estimation (internal): starting values and screening |
| `R/cnorm.R`, `R/cgnorm.R` | complex normal; complex generalised normal (`dcgnorm`, `rcgnorm`, internal `cgnormConcentrated`) |
| `R/cacf.R`, `R/cvar.R` | complex ACF/PACF, direct/conjugate variances and correlations |
| `R/cTransformations.R` | `cscale`, `cdescale`, `clog`, `cexp` |
| `src/complexCode.cpp` | Rcpp: `invert`, `polyprodcomplex` |
| `tests/testthat/` | tests (none existed before v1.0.3) |
| `dev/` | plans and prototypes (in `.Rbuildignore`); `dev/PLAN-dcgnorm.md` has the design decisions |

## Workflow

```bash
cd ~/R/Projects/Packages/complex
Rscript -e 'roxygen2::roxygenise()'          # roxygen2 8.1.0; never edit NAMESPACE/man by hand
R CMD INSTALL .                               # installs into the user library
Rscript -e 'library(complex); testthat::test_dir("tests/testthat", reporter="summary")'
R CMD build --no-build-vignettes . && R CMD check --no-manual --ignore-vignettes complex_*.tar.gz
```
- Build and check in a scratch directory, not in the package folder.
- Expected check result: `Status: 1 WARNING`, the "OS reports request to set locale to en_US.UTF-8"
  warning of this machine. Anything else must be fixed.
- Add a NEWS entry (Changes / Bugfixes) for every user-visible change.
- `src/*.o` and `src/complex.so` are tracked in git (historical). Do not commit changes to them.

## Code style

Match the existing code: 4-space indentation, `;` at the end of statements, `if(...){ ... }` with
`else{` on the next line after `}`, camelCase names, `<-`, comments describing why. Use
`seq_len(n)` rather than `1:n`: several bugs here came from `1:0` selecting the first element.

## clm() internals worth knowing

- Parameters `B` are complex (intercept, regressors, AR, MA). The optimiser (nloptr, SBPLX by default)
  works on `c(Re(B), Im(B))`; for `distribution="dcgnorm"` with estimated shape, log(shape) is
  appended at the end. A starting shape can be passed as the last element of `B`.
- The ARI part is stored via the polynomial `(1-phi(B))(1-B)^d`; `fitter()` builds the full
  parameter vector. Initial lags are extrapolated (`xregExpander(gaps="auto")`), initial MA errors
  are zero; all observations enter the likelihood.
- Likelihood `dcnorm`: concentrated, `-T(log 2pi + 1 + 0.5 log det Sigma)`. `dcgnorm`: circularity
  from moments, scale by exact ML given it (approximate ML for shape != 1).
- The loss returns `1E+100/min|root|` if the complex AR or MA polynomial has a root inside the unit
  circle (differencing is not checked).
- `nparam()` counts per series (complex units: covariance 3/2, shape 1/2). `logLik` df and the ICs
  use all real parameters (2*nparam); `clmIC()` has the formulas (AICc/BICc: Bedrick & Tsai with
  two series, approximate for the restricted complex model).
- `fast=TRUE`: skips data checks and, for pure cARIMA, starts from Hannan-Rissanen estimates.
- `vcov.clm` re-estimates the model from `object$call` with `parameters` fixed and `FI=TRUE`; the
  MA columns of `object$data` then act as regressors, so the parameter vector does not match the
  orders there. Keep `object$call` a clm call (carima stores its own call in `carimaCall`).
- `loss="CLS"` cannot be used with MA (complex-valued loss).

## auto.carima()

- stepwise (default): for each d, neighbourhood search over p, q (IC only); constant included for
  d=0, excluded for d>0, then flipped for the best model of each d and the search continued.
- Every candidate starts from Hannan-Rissanen and from the closest estimated model.
- ICs are computed on the common sample without the first `ar+i+ma` observations (models are
  estimated on the full sample).
- `fast=TRUE`: Hannan-Rissanen screening of all models, likelihood refit of the best per
  (d, constant, q), then neighbourhood search. About 3x faster, but approximate (least squares
  underrates MA models with heavy tails).
- `search="full"`: all combinations.

## Open items

- Phase 3/4/5 of `dev/PLAN-dcgnorm.md`: joint Hessian with the shape, simulated prediction
  intervals, dcgnorm diagnostics.
- Seasonal cARIMA is not supported.
- A complex portmanteau / partial autocorrelation statistic combining direct and conjugate ACFs
  (discussed, not implemented).
