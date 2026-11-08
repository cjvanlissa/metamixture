# Calculate Surprisal Using Kernel Density Estimation

Calculates the surprisal (negative log density) of observed values under
a distribution estimated from a reference sample using kernel density
estimation (KDE).

## Usage

``` r
surprisal(x, y, aggregate = "sum", ...)
```

## Arguments

- x:

  Numeric vector containing the reference sample used to estimate the
  probability density.

- y:

  Numeric vector containing the values at which surprisal is evaluated.

- aggregate:

  Either `"sum"` (default), `"mean"`, or `FALSE`. Determines whether
  surprisal values are summed, averaged, or returned individually.

- ...:

  Additional arguments passed to
  [`stats::density()`](https://rdrr.io/r/stats/density.html).

## Value

A numeric vector containing individual surprisal values if
`aggregate = FALSE`, or a single numeric value containing their sum or
mean otherwise.

## Examples

``` r
set.seed(123)

x <- rnorm(1000)
y <- c(0, 1, 2, 3)

# Individual surprisal values
surprisal(x, y, aggregate = FALSE)
#> [1] 0.8767707 1.4248554 2.7872109 5.6401511

# Total surprisal
surprisal(x, y)
#> [1] 10.72899

# Mean surprisal
surprisal(x, y, aggregate = "mean")
#> [1] 2.682247
```
