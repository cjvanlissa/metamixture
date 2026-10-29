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

#' Posterior class probabilities from component log-densities and weights
#' (works for any mixture): logdens is n x K, weights is length K.
posterior_from_logdens <- function(logdens, weights) {
  lt <- sweep(as.matrix(logdens), 2, log(weights), "+")
  m <- apply(lt, 1, max)
  lse <- m + log(rowSums(exp(lt - m)))
  exp(lt - lse)
}



icl_mclust <- function(fit) icl(fit$z, fit$loglik, fit$df, fit$n)

icl_mx <- function(model) {
  z <- class_prob_individual_remix(model)
  ll <- -0.5 * as.numeric(model$output$fit)
  npar <- length(model$output$estimate)   # counts free parameters (respects equality labels)
  icl(z, ll, npar, n = length(model$Y$values))
}

class_prob_individual_remix <- function(model){
  yi <- as.numeric(model$Y$values)
  vi <- as.numeric(model$V$values)
  mu <- model@matrices$M$values[1, , drop = TRUE]
  tau2 <- diag(model@matrices$S$values)
  K <- length(mu)
  w <- if (K == 1) 1 else as.numeric(OpenMx::mxEval(weights, model))
  logd <- sapply(seq_len(K), function(k) stats::dnorm(yi, mu[k], sqrt(vi + tau2[k]), log = TRUE))
  logd <- matrix(logd, ncol = K)
  posterior_from_logdens(logd, w)
}
