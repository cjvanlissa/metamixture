# Cross-Validated Class Enumeration

Selects the number of latent classes in random-effects mixture
meta-analysis using K-fold cross-validation. Candidate models are
compared based on their out-of-sample negative log predictive density
(NLPD).

## Usage

``` r
cv_enumeration(x, folds = 5L)
```

## Arguments

- x:

  An object for which a method exists, typically a `mixture_list` of at
  least two fitted `remix` models with different numbers of latent
  classes. The models must be based on the same observations and use the
  same constraints on class-specific means and heterogeneity variances.
  Models should be ordered by increasing number of classes, starting
  with a one-class model.

- folds:

  Integer. Number of cross-validation folds. Defaults to 5. Must be
  between 2 and the number of observations.

## Value

A list containing:

- best_class:

  Integer index of the candidate model with the lowest total
  cross-validated negative log predictive density. If models are ordered
  from one to the maximum number of classes, this index corresponds to
  the selected number of classes.

- scores:

  Named numeric vector containing the total negative log predictive
  density across folds for each candidate model.

- fold_scores:

  Numeric matrix containing the negative log predictive density for each
  fold (rows) and candidate model (columns).

- fold_id:

  Integer vector indicating the fold assigned to each observation.

## Details

The function randomly partitions the observations into `folds`
approximately equal-sized subsets. For each candidate number of classes,
the model is re-estimated on the training observations and evaluated on
the held-out observations using the negative log predictive density.
Smaller values indicate better out-of-sample predictive performance.

The negative log predictive densities are summed across folds for each
candidate model. The model with the smallest total score is selected.

The mean and heterogeneity variance constraints are inferred from the
second model in `x` and applied consistently to all models fitted during
cross-validation.

Because fold assignments are random, results may differ between calls.
Use [`set.seed()`](https://rdrr.io/r/base/Random.html) before calling
the function to obtain reproducible fold assignments and random starting
values.

## Examples

``` r
if (FALSE) { # \dontrun{
library(metamixture)

set.seed(123)

# Simulate meta-analytic data
dat <- sim_remix_fixed_tau(n = 100, es = 2)

# Estimate candidate models
models <- lapply(1:3, function(k) {
remix(
  yi = dat$yi,
  vi = dat$vi,
  classes = k,
  tau2 = "equal"
)
})

# Five-fold cross-validation
cv <- cv_enumeration(models, folds = 5)

# Selected number of classes
cv$best_class

# Total predictive scores
cv$scores

# Scores by fold and candidate model
cv$fold_scores
} # }
```
