# Random-Effects Mixture Meta-Analysis

Fits finite mixture models for random-effects meta-analysis using
maximum likelihood estimation in `OpenMx`. Each mixture class has its
own mean effect size and residual heterogeneity variance. Note that
`OpenMx` must be loaded using
[`library(OpenMx)`](https://openmx.ssri.psu.edu/) before calling this
function.

## Usage

``` r
remix(
  yi,
  vi,
  classes = 2L,
  mu = c("varying", "equal"),
  tau2 = c("varying", "equal"),
  run = TRUE,
  ...
)
```

## Arguments

- yi:

  Numeric vector of observed effect sizes.

- vi:

  Numeric vector of known sampling variances.

- classes:

  Integer vector specifying the number of mixture classes to fit.
  Defaults to `1L:2L`.

- mu:

  Whether component means are `varying` or `equal`.

- tau2:

  Whether component heterogeneity variances are `varying` or `equal`.

- run:

  Logical. If `TRUE` (the default), the model is estimated using OpenMx.
  If `FALSE`, the unestimated OpenMx model object is returned, allowing
  further modification before estimation.

- ...:

  Arguments passed to other functions.

## Value

A fitted model object, or a `"mixture_list"` when multiple class counts
are specified.

## Details

The estimation procedure is described in Van Lissa, Garnier-Villarreal,
& Anadria (2023): starting values are determined by way of K-means
clustering, then the model is estimated using simulated annealing,
followed by maximum likelihood optimization.

## References

Van Lissa, C. J., Garnier-Villarreal, M., & Anadria, D. (2023).
Recommended Practices in Latent Class Analysis using the Open-Source
R-Package tidySEM. Structural Equation Modeling.
[doi:10.1080/10705511.2023.2250920](https://doi.org/10.1080/10705511.2023.2250920)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(7894)
df <- sim_remix_fixed_tau()

# Conventional random-effects meta-analysis
res1 <- remix(df$yi, df$vi, classes = 1)
if(requireNamespace("metafor")){
  res_mf <- metafor::rma(yi = df$yi, vi = df$vi)
  res_mf$b[1,1] - res1$M$values[1, 1]
  res_mf$tau2 - res1$S$values[1, 1]
}

# Two-class random-effects mixture model
res2 <- remix(df$yi, df$vi, classes = 2)
table_fit(res2)
table_results(res2)
} # }
```
