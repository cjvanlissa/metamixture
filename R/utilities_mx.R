remix_weights <- function(x){
  if(ncol(x$M$values) == 1L){
    return(1)
  }
  x$weights$result[1, , drop = TRUE]
}

num_classes_remix <- function(x){
  ncol(x@matrices$M$values)
}

# Calculate probabilities from assigned model parameters
p_from_parameters <- function(x){
  yi <- x$Y$values
  vi <- x$V$values
  classes <- num_classes(x)
  mu <- x@matrices$M$values[1,]
  tau2 <- diag(x@matrices$S$values)

  # Log of weighted component densities
  logdens <- vapply(seq_len(classes), function(k) {
    stats::dnorm(yi, mean = mu[k], sd = sqrt(vi + tau2[k]), log = TRUE)
  }, numeric(length(yi)))
  logdens <- matrix(logdens, nrow = length(yi), ncol = classes)
  row_max <- apply(logdens, 1, max)
  posteriors <- exp(logdens - row_max)
  return(colMeans(posteriors/rowSums(posteriors)))
}
