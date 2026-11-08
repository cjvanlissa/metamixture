# Class Enumeration Using Predictive Surprisal

Compares random-effects mixture meta-analysis models using
simulation-based predictive surprisal. For each candidate model,
repeated datasets are simulated, and kernel density estimation (KDE) is
used to evaluate the surprisal of the observed effect sizes under the
simulated distribution.

## Usage

``` r
pmc_surprisal(x, ..., reps = 100L, ci = 0.95, nboot = 100L)
```

## Arguments

- x:

  A list of at least two fitted `remix` models, typically representing
  different numbers of latent classes. All models must be fitted to the
  same observed effect sizes.

- ...:

  Additional arguments passed to
  [`surprisal()`](https://cjvanlissa.github.io/metamixture/reference/surprisal.md).

- reps:

  Integer. Number of simulated datasets generated for each candidate
  model. Defaults to 100.

- ci:

  Numeric. Confidence level for the bootstrap confidence intervals of
  differences in mean surprisal. Defaults to 0.95.

- nboot:

  Integer. Number of bootstrap replications used to estimate confidence
  intervals for differences in mean surprisal. Defaults to 100.

## Value

A data.frame of class `pmc_df` containing:

- comparison:

  Comparison type: `dif_seq` or `dif_one`.

- null:

  Name of the reference model.

- alt:

  Name of the alternative model.

- null_stat:

  Mean simulated-data surprisal for the reference model.

- alt_stat:

  Mean simulated-data surprisal for the alternative model.

- difference:

  Difference in mean surprisal (`alt_stat - null_stat`).

- lb:

  Lower bootstrap confidence limit.

- ub:

  Upper bootstrap confidence limit.

- sig:

  An asterisk (`*`) if the upper confidence limit is below zero;
  otherwise an empty string.

## Details

For each candidate model, the function generates `reps` simulated
datasets.

For each simulated dataset,
[`surprisal()`](https://cjvanlissa.github.io/metamixture/reference/surprisal.md)
estimates the distribution of the simulated effect sizes using kernel
density estimation and evaluates the surprisal of the observed effect
sizes under this distribution.

The mean surprisal across simulated datasets is used as a measure of
predictive discrepancy. Smaller values indicate that the simulated
distributions assign greater density to the observed effect sizes and
therefore provide a better approximation to the observed data.

Two types of model comparisons are performed:

- dif_seq:

  Compares each model with the preceding model in the supplied list.

- dif_one:

  Compares each model with the first model in the supplied list.

For each comparison, the difference is calculated as the mean surprisal
of the alternative model minus the mean surprisal of the reference
model:

\$\$\Delta = \bar{S}\_{alt} - \bar{S}\_{null}\$\$

Negative differences favor the alternative model.

Confidence intervals are obtained using a nonparametric percentile
bootstrap. The simulated surprisal scores for the two models are
independently resampled with replacement, and the difference in their
means is recalculated for each bootstrap replication.

A comparison is marked with an asterisk (`*`) when the upper confidence
limit is below zero, indicating that the alternative model has
significantly lower mean surprisal at the specified confidence level.

If the `future.apply` package is installed, simulations are evaluated
using
[`future.apply::future_replicate()`](https://future.apply.futureverse.org/reference/future_lapply.html),
allowing parallel execution when an appropriate `future` plan has been
configured. Otherwise, the function uses sequential replication.

The function uses random simulation and bootstrap resampling. Use
[`set.seed()`](https://rdrr.io/r/base/Random.html) for reproducibility.

## Examples

``` r
if (FALSE) { # \dontrun{
library(metamixture)

set.seed(123)

# Simulate meta-analytic data
dat <- sim_remix_fixed_tau(k = 100, es = 2)

# Fit candidate models
models <- lapply(1:3, function(k) {
remix(
  yi = dat$yi,
  vi = dat$vi,
  classes = k,
  tau2 = "equal"
)
})

# Compare models using predictive surprisal
pmc <- pmc_surprisal(
models,
reps = 100,
nboot = 500
)

# Model comparisons
pmc

# Simulated surprisal scores
head(attr(pmc, "scores"))
} # }

```
