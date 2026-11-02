# Posterior class probabilities from component log-densities and weights
posterior_from_logdens <- function(logdens, weights) {
  lt <- sweep(as.matrix(logdens), 2, log(weights), "+")
  m <- apply(lt, 1, max)
  lse <- m + log(rowSums(exp(lt - m)))
  exp(lt - lse)
}
extract_postprob_remix <- function(model){
  yi <- as.numeric(model$Y$values)
  vi <- as.numeric(model$V$values)
  mu <- model@matrices$M$values[1, , drop = TRUE]
  tau2 <- diag(model@matrices$S$values)
  K <- length(mu)
  w <- if (K == 1) 1 else as.numeric(OpenMx::mxEvalByName("weights", model))
  logd <- sapply(seq_len(K), function(k) stats::dnorm(yi, mu[k], sqrt(vi + tau2[k]), log = TRUE))
  logd <- matrix(logd, ncol = K)
  posterior_from_logdens(logd, w)
}
