surprisal <- function(x, y, ...) {
  # Estimate density of x using KDE
  kde <- density(x, ...)

  # Evaluate estimated density at each y
  # rule = 2 extends boundary values outside KDE grid
  p <- approx(kde$x, kde$y, xout = y, rule = 2)$y

  # Calculate surprisal
  # Negative log kernel density
  return(-log(pmax(p, .Machine$double.xmin)))
}
