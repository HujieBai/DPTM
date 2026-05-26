# DPTM

`DPTM` is an R package for estimating dynamic panel threshold models with fixed
effects. It supports both dynamic panel linear models and dynamic panel
multiple-threshold models, with tools for threshold estimation, model
comparison, bootstrap testing, and convergence diagnostics.

The package is designed for applied panel-data work where the effect of one or
more regressors may change across regimes determined by an observed threshold
variable.

## Features

- Estimate fixed-effects dynamic panel linear models with `DPML()`.
- Estimate dynamic panel threshold models with `DPTS()`.
- Allow one or more thresholds through the `Th` argument.
- Separate regressors with threshold-varying effects from controls with common
  effects.
- Include or exclude threshold effects for the lagged dependent variable.
- Add time fixed effects when needed.
- Search unknown thresholds by MCMC-MLE or deterministic grid search.
- Use joint grid search, sequential grid search, and coordinate refinement for
  multiple-threshold specifications.
- Test threshold existence and the number of thresholds with bootstrap
  likelihood-ratio-type statistics.
- Report coefficient estimates, standard errors, z statistics, threshold
  estimates, likelihood values, and MCMC diagnostics.

## Installation

Install from GitHub after replacing `<owner>` with the account or organization
that hosts the repository:

```r
install.packages("remotes")
remotes::install_github("<owner>/DPTM")
```

You can also install a local source archive:

```r
install.packages("DPTM_4.0.1.tar.gz", repos = NULL, type = "source")
```

Or install from a local package directory:

```r
remotes::install_local("path/to/DPTM")
```

## Quick Start

```r
library(DPTM)

data(d1)
head(d1)
```

The example dataset `d1` is a simulated balanced panel with individual and time
indexes, a dependent variable, a lagged dependent variable, a threshold variable,
a regressor with threshold effects, and a control variable.

### Dynamic Panel Linear Model

Use `DPML()` to estimate the linear dynamic panel model with fixed effects.

```r
linear_model <- DPML(
  y ~ x + z,
  data = d1,
  index = c("id", "year")
)

print(linear_model)
```

Add time fixed effects with `timeFE = TRUE`.

```r
linear_time_fe <- DPML(
  y ~ x + z,
  data = d1,
  index = c("id", "year"),
  timeFE = TRUE
)
```

### One-Threshold Model

Use `DPTS()` to estimate a dynamic panel threshold model. The first formula
contains regressors with threshold-varying effects. The second formula contains
controls whose effects are common across regimes.

```r
set.seed(42)

one_threshold <- DPTS(
  y ~ x,
  y ~ z,
  data = d1,
  index = c("id", "year"),
  q = d1$q,
  Th = 1,
  iterations = 1000
)

print(one_threshold)
```

The printed result includes the threshold estimate, regime-specific
coefficients, the negative log-likelihood value, confidence intervals for
thresholds when MCMC search is used, and Gelman-Rubin diagnostics.

### Multiple-Threshold Model

Set `Th > 1` to estimate multiple thresholds.

```r
set.seed(42)

two_thresholds <- DPTS(
  y ~ x,
  y ~ z,
  data = d1,
  index = c("id", "year"),
  q = d1$q,
  Th = 2,
  iterations = 1000
)

print(two_thresholds)
```

This specification partitions the sample into three regimes:

- `q <= gamma_1`
- `gamma_1 < q <= gamma_2`
- `q > gamma_2`

## Choosing the Threshold Search Method

By default, `DPTS()` searches thresholds with an MCMC-MLE procedure based on
`BayesianTools`. This is useful when the number of thresholds is greater than
one, since exhaustive grid search can become expensive.

```r
model_mcmc <- DPTS(
  y ~ x,
  y ~ z,
  data = d1,
  index = c("id", "year"),
  q = d1$q,
  Th = 2,
  iterations = 2000
)
```

For deterministic threshold search, set `grid_search = TRUE`.

```r
model_grid <- DPTS(
  y ~ x,
  y ~ z,
  data = d1,
  index = c("id", "year"),
  q = d1$q,
  Th = 1,
  grid_search = TRUE,
  grids = 100
)
```

For multiple thresholds, choose between joint and sequential grid search.

```r
model_joint <- DPTS(
  y ~ x,
  y ~ z,
  data = d1,
  index = c("id", "year"),
  q = d1$q,
  Th = 2,
  grid_search = TRUE,
  grids = 100,
  grid_search_type = "jointly",
  grid_search_iter = 1
)

model_sequential <- DPTS(
  y ~ x,
  y ~ z,
  data = d1,
  index = c("id", "year"),
  q = d1$q,
  Th = 2,
  grid_search = TRUE,
  grids = 100,
  grid_search_type = "sequential",
  grid_search_iter = 1
)
```

`grid_search_iter` controls coordinate-refinement cycles after the initial
search. For two thresholds, one refinement cycle re-estimates the first
threshold conditional on the second and then re-estimates the second conditional
on the first.

## Testing Threshold Effects

`Threshold_Test()` implements bootstrap tests for threshold existence and for
the number of thresholds.

```r
set.seed(42)

test_one <- Threshold_Test(
  y ~ x,
  y ~ z,
  data = d1,
  index = c("id", "year"),
  q = d1$q,
  Th = 1,
  bt = 50,
  iterations = 500
)

test_one
```

Interpretation of `Th` in `Threshold_Test()`:

- `Th = 1`: test the null of no threshold against one threshold.
- `Th = 2`: test the null of one threshold against two thresholds.
- `Th = 3`: test the null of two thresholds against three thresholds.

The same grid-search options used by `DPTS()` are available in
`Threshold_Test()`.

```r
test_grid <- Threshold_Test(
  y ~ x,
  y ~ z,
  data = d1,
  index = c("id", "year"),
  q = d1$q,
  Th = 1,
  bt = 50,
  grid_search = TRUE,
  grids = 100,
  grid_search_type = "sequential"
)
```

## Model Specification Guide

`DPTS()` uses two formula arguments:

- `formula`: regressors that may have regime-specific coefficients.
- `formula_cv`: regressors that enter as controls with common coefficients.

Examples:

```r
# Threshold effects in x; common effect for z
DPTS(y ~ x, y ~ z, data = d1, index = c("id", "year"), q = d1$q)

# Threshold effects in x, but no threshold effects in lagged y
DPTS(y ~ x, y ~ z, data = d1, index = c("id", "year"), q = d1$q, NoY = TRUE)

# Threshold effects only in lagged y; x and z are controls
DPTS(NULL, y ~ x + z, data = d1, index = c("id", "year"), q = d1$q)
```

Useful arguments:

- `index`: character vector naming the individual and time indexes.
- `q`: threshold variable.
- `Th`: number of thresholds in the model.
- `timeFE`: whether to include time fixed effects.
- `NoY`: whether the lagged dependent variable should enter without threshold
  effects.
- `y1`: user-supplied lagged dependent variable, useful for transformed
  dependent variables.
- `sro`: minimum regime proportion used to restrict threshold candidates.
- `r0x`, `r1x`: lower and upper bounds of the threshold search interval.
- `iterations`: number of MCMC iterations; half are treated as burn-in.
- `grid_search`, `grids`, `grid_search_type`, `grid_search_iter`: deterministic
  search controls.

## Returned Objects

`DPML()` and `DPTS()` return objects of class `"DPTM"`. Important fields include:

- `coefficients`: estimated coefficients.
- `NNLL`: negative log-likelihood value.
- `Zvalues`: z statistics.
- `Ses`: standard errors.
- `covariance_matrix`: covariance matrix of estimated coefficients.
- `Th`: number of thresholds.
- `thresholds`: estimated thresholds.
- `threshold_search`: `"MCMC"` or `"grid"`.
- `grid_results`: searched threshold candidates and objective values when grid
  search is used.

`Threshold_Test()` returns an `"htest"` object with:

- `statistic`: bootstrap test statistic.
- `p.value`: bootstrap p-value.
- `estimate`: critical value.
- `LRs`: bootstrap statistics.

## Included Data

- `d1`: simulated dynamic panel data for package examples.
- `Growth_Inflation`: country-level growth and inflation panel data with
  five-year periods.

## Practical Notes

- MCMC search can be computationally intensive. Increase `iterations` for final
  analysis and use smaller values only for testing code.
- Grid search is deterministic and easier to reproduce, but can become expensive
  for multiple thresholds.
- For multiple thresholds, `grid_search_type = "sequential"` is often faster
  than `grid_search_type = "jointly"`.
- MCMC threshold confidence intervals and Gelman-Rubin diagnostics are available
  only when MCMC search is used.
- Bootstrap tests can be slow. Increase `bt` for final inference.

## License

This package is released under the GPL (>= 3) license.
