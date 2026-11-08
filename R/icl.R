#' Integrated Completed Likelihood (ICL) for any mixture model
#'
#' ICL = BIC + 2 * Entropy, on the "lower is better" (-2LL) scale, where
#' Entropy = -sum_ik z_ik * log(z_ik) (Biernacki, Celeux & Govaert, 2000).
#'
#' @param z  n x K matrix of posterior class probabilities (rows sum to 1).
#' @param loglik Log-likelihood of the fitted model.
#' @param npar Number of free parameters.
#' @param n Number of observations (default: nrow(z)).
#' @return Named numeric vector: loglik, npar, n, BIC, entropy, ICL.
icl <- function(z, loglik, npar, n = nrow(z)) {
  tryCatch({
    z <- as.matrix(z)
  if (any(z < -1e-8, na.rm = TRUE) || any(abs(rowSums(z) - 1) > 1e-6, na.rm = TRUE)) {
    message("`z` must be a matrix of posterior probabilities with rows summing to 1.")
  }
  zlogz <- ifelse(z > 0, z * log(z), 0)   # define 0 * log(0) = 0
  entropy <- -sum(zlogz)
  bic <- -2 * loglik + npar * log(n)
  bic + 2 * entropy}, error = function(e){NA})
}
