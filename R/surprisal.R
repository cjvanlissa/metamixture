#' Calculate Surprisal Using Kernel Density Estimation
#'
#' Calculates the surprisal (negative log density) of observed values
#' under a distribution estimated from a reference sample using
#' kernel density estimation (KDE).
#'
#' @param x Numeric vector containing the reference sample used to estimate the
#' probability density.
#' @param y Numeric vector containing the values at which surprisal is
#' evaluated.
#' @param aggregate Either `"sum"` (default), `"mean"`, or `FALSE`.
#' Determines whether surprisal values are summed, averaged,
#' or returned individually.
#' @param ... Additional arguments passed to [stats::density()].
#' @return
#' A numeric vector containing individual surprisal values
#' if `aggregate = FALSE`, or a single numeric value
#' containing their sum or mean otherwise.
#'
#' @examples
#' set.seed(123)
#'
#' x <- rnorm(1000)
#' y <- c(0, 1, 2, 3)
#'
#' # Individual surprisal values
#' surprisal(x, y, aggregate = FALSE)
#'
#' # Total surprisal
#' surprisal(x, y)
#'
#' # Mean surprisal
#' surprisal(x, y, aggregate = "mean")
#' @export
surprisal <- function(x, y, aggregate = "sum", ...) {
  # Estimate density of x using KDE
  kde <- stats::density(x, ...)

  # Evaluate estimated density at each y
  # rule = 2 extends boundary values outside KDE grid
  p <- stats::approx(kde$x, kde$y, xout = y, rule = 2)$y

  # Calculate surprisal
  # Negative log kernel density
  out <- -log(pmax(p, .Machine$double.xmin))
  if(isFALSE(aggregate)){
    return(out)
  }
  if(identical(aggregate, "mean")) return(mean(out))
  return(sum(out))
}
