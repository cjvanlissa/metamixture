# Cluster-Robust Standard Errors for `remix` Models

Computes cluster-robust sandwich standard errors and confidence
intervals for random-effects mixture models fitted using
[`remix()`](https://cjvanlissa.github.io/metamixture/reference/remix.md).
This allows for dependence between effect sizes belonging to the same
cluster, such as multiple effect sizes from the same study.

## Usage

``` r
sandwich_remix(model, cluster, type = c("CR1", "CR0"), rel_step = 1e-04)
```

## Arguments

- model:

  A fitted OpenMx model created using
  [`remix()`](https://cjvanlissa.github.io/metamixture/reference/remix.md).

- cluster:

  Vector identifying the cluster membership of each effect size, e.g.,
  the study the effect size came from. Must have the same length and
  order as the observed effect sizes.

- type:

  Cluster-robust covariance estimator. Either `"CR1"` (default),
  applying the correction G/(G-1), or `"CR0"`, the unadjusted sandwich
  estimator, where G is the number of clusters (or, studies).

- rel_step:

  Relative step size for numerical differentiation of the analytic score
  function. Defaults to `1e-4`.

## Value

An `MxModel` with `attr("se_robust")` containing the standard errors,
and `attr("sandwich")` with `vcov`, the cluster-robust covariance
matrix, `vcov_free`, the cluster-robust covariance matrix for the unique
free model parameters, `information`, observed Fisher information matrix
(negative Hessian of the log-likelihood), `cluster_scores`, the matrix
of first derivatives of the log-likelihood for each cluster (rows) by
free parameters (columns); `n` number of observed effect sizes,
`clusters` number clusters G, `npar`, number of freely estimated
parameters, `type` passed from input.

## Details

Note that cluster-robust standard errors provide consistent inference
for dependent effect sizes but do not affect the fitted parameter
estimates, nor model the dependency in the data. CR0 and CR1 are
large-sample approximations and may be unreliable when the number of
clusters is small or mixture parameters are weakly identified.

## Examples

``` r
if (FALSE) { # \dontrun{
library(metamixture)

set.seed(1)

n <- 100
df <- metamixture::sim_remix_fixed_tau(n, es = 2)

res1 <- remix(
  yi = df$yi,
  vi = df$vi,
  classes = 2,
  tau2 = "equal"
)

# Original sample: one observation per cluster
sw1 <- sandwich_remix(
  res1,
  cluster = seq_len(n),
  type = "CR0"
)

# Duplicate observations exactly
df2 <- rbind(df, df)
df2$study <- rep(seq_len(n), 2)

res2 <- remix(
  yi = df2$yi,
  vi = df2$vi,
  classes = 2,
  tau2 = "equal"
)

sw2 <- sandwich_remix(
  res2,
  cluster = df2$study,
  type = "CR0"
)

# Compare robust covariance matrices
all.equal(
  attr(sw1, "sandwich")$vcov,
  attr(sw2, "sandwich")$vcov,
  tolerance = 1e-3
)
} # }
```
