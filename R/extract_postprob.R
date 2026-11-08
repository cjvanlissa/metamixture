# Convert component log-densities to posterior probabilities
posterior_from_logdens <- function(logdens, weights) {

  logdens <- as.matrix(logdens)
  weights <- as.numeric(weights)

  # Log of weighted component densities
  lt <- sweep(logdens, 2, log(weights), "+")

  # Stable log-sum-exp normalization
  m <- apply(lt, 1, max)

  centered <- sweep(lt, 1, m, "-")
  denom <- rowSums(exp(centered))

  exp(sweep(centered, 1, log(denom), "-"))
}


# Class-specific densities or posterior probabilities
class_densities <- function(model, yi, vi, normalize = TRUE) {
  mu <- as.numeric(OpenMx::mxEvalByName("M", model))
  tau2 <- diag(OpenMx::mxEvalByName("S", model))
  K <- length(mu)

  weights <- if (K == 1L) {
    1
  } else {
    as.numeric(OpenMx::mxEvalByName("weights", model))
  }

  # Component log-densities, without mixing weights
  logdens <- vapply(seq_len(K), function(k) {

    stats::dnorm(
      yi,
      mean = mu[k],
      sd = sqrt(vi + tau2[k]),
      log = TRUE
    )

  }, numeric(length(yi)))

  logdens <- matrix(logdens, ncol = K)

  colnames(logdens) <- paste0("class", seq_len(K))

  if (!normalize) {
    return(logdens)
  }

  posterior_from_logdens(logdens, weights)
}

extract_postprob_remix <- function(model) {

  class_densities(
    model = model,
    yi = as.numeric(model$Y$values),
    vi = as.numeric(model$V$values),
    normalize = TRUE
  )
}
