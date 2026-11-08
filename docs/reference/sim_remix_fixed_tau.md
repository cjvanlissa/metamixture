# Simulate Data from a Random-Effects Mixture Model

Simulates effect sizes from a finite mixture of normal random-effects
models, with class-specific means and residual heterogeneity.

## Usage

``` r
sim_remix_fixed_tau(
  k = 100L,
  classes = 2L,
  p_ratio = 1,
  vi_range = c(0.001, 0.01),
  es = 1
)
```

## Arguments

- k:

  Number of effect sizes to simulate.

- classes:

  Number of mixture classes.

- p_ratio:

  Numeric, balance of class proportions. Default `1` equally balances
  class proportions; higher values assign a higher proportion to the
  first class; e.g., `2` means the first class will have twice the
  proportion of cases as other classes.

- vi_range:

  Numeric vector specifying the minimum and maximum sampling variances.

- es:

  Numeric, effect size (distance between class centroids).

## Value

A data frame containing the simulated effect sizes (`yi`), sampling
variances (`vi`), and true class memberships (`C`).

## Examples

``` r
set.seed(123)
dat <- sim_remix_fixed_tau(k = 5)
dat
#>           yi   vi C
#> 1 -0.6979824 0.01 2
#> 2  1.2456779 0.01 1
#> 3  0.8904906 0.01 2
#> 4 -0.1178267 0.01 1
#> 5  0.1839957 0.01 1
```
