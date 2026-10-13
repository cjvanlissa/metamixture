surprisal <- function(x, y, aggregate = "sum", ...) {
  # Estimate density of x using KDE
  kde <- density(x, ...)

  # Evaluate estimated density at each y
  # rule = 2 extends boundary values outside KDE grid
  p <- approx(kde$x, kde$y, xout = y, rule = 2)$y

  # Calculate surprisal
  # Negative log kernel density
  out <- -log(pmax(p, .Machine$double.xmin))
  if(isFALSE(aggregate)){
    return(out)
  }
  if(identical(aggregate, "mean")) return(mean(out))
  return(sum(out))
}
