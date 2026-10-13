#' Observed-information standard errors for a remix fit
#'
#' Computes the inverse observed information using finite differences of the
#' analytic observed-data score. Requires an interior, identifiable ML solution.
#' @param fit Result from memi::remix().
#' @param yi Observed effect sizes used to fit the model.
#' @param vi Known sampling variances used to fit the model.
#' @param rel_step Relative finite-difference step on scaled parameters.
#' @param boundary_tol Tolerance for recognizing boundary estimates.
#' @return List with se_mu, se_tau2, se_pi, vcov, information, estimates,
#'   and the minimum eigenvalue of the scaled information matrix.
#' @export
remix_se <- function(fit, yi, vi, rel_step = 1e-4,
                     boundary_tol = 1e-10) {
  yi <- as.numeric(yi); vi <- as.numeric(vi)
  K <- length(fit$mu)
  if (length(yi) != length(vi) || !length(yi) ||
      any(!is.finite(yi)) || any(!is.finite(vi)) || any(vi <= 0))
    stop("yi and vi must be equally long, finite, with vi > 0.")
  if (length(fit$tau2) != K || length(fit$pi) != K ||
      any(!is.finite(c(fit$mu, fit$tau2, fit$pi))))
    stop("Invalid fitted parameters.")
  if (!isTRUE(fit$converged))
    warning("Fit has not reported convergence; standard errors may be unreliable.")
  if (any(fit$tau2 <= boundary_tol * max(1, stats::var(yi))) ||
      any(fit$pi <= boundary_tol) || (K > 1L && any(fit$pi >= 1 - boundary_tol)))
    stop("Boundary fit: inverse-information SEs are not reliable. Use a parametric bootstrap or profile likelihood.")
  if (K > 1L && any(diff(sort(fit$mu)) == 0))
    warning("Coincident component means; inspect identifiability.")
  if (!is.finite(rel_step) || rel_step <= 0 || rel_step >= 0.1)
    stop("rel_step must be in (0, 0.1).")
  n <- length(yi)
  # Centering and scaling improve the numerical differentiation.
  center <- mean(yi)
  scale <- max(stats::sd(yi), sqrt(stats::median(vi)), 1e-8)
  y <- (yi - center) / scale
  v <- vi / scale^2
  m <- as.numeric((fit$mu - center) / scale)
  t <- as.numeric(fit$tau2 / scale^2)
  pi <- as.numeric(fit$pi)
  eta <- if (K > 1L) log(pi[seq_len(K - 1L)] / pi[K]) else numeric()
  theta <- c(m, t, eta)
  p <- length(theta)
  softmax <- function(eta) {
    z <- c(eta, 0)
    z <- exp(z - max(z))
    z / sum(z)
  }
  score <- function(th) {
    mu <- th[seq_len(K)]
    tau <- th[K + seq_len(K)]
    probs <- softmax(if (K > 1L) th[2L * K + seq_len(K - 1L)] else numeric())
    if (any(tau <= -min(v))) stop("Nonpositive component variance in numerical differentiation.")
    s <- outer(v, tau, `+`)
    resid <- outer(y, mu, `-`)
    logd <- sweep(-0.5 * (log(2 * base::pi) + log(s) + resid^2 / s),
                  2L, log(probs), `+`)
    rowmax <- apply(logd, 1L, max)
    a <- exp(logd - rowmax)
    r <- a / rowSums(a)
    smu <- colSums(r * resid / s)
    stau <- colSums(0.5 * r * (resid^2 / s^2 - 1 / s))
    seta <- if (K > 1L) colSums(r)[seq_len(K - 1L)] - n * probs[seq_len(K - 1L)] else numeric()
    c(smu, stau, seta)
  }
  # Differentiate score with respect to each parameter. Use an adaptive
  # centered step, keeping tau2 perturbations strictly inside the domain.
  H <- matrix(0, p, p)
  for (j in seq_len(p)) {
    h <- rel_step * max(1, abs(theta[j]))
    if (j > K && j <= 2L * K) h <- min(h, theta[j] * 0.25)
    if (h <= 0 || !is.finite(h)) stop("Could not construct a valid finite-difference step.")
    plus <- minus <- theta
    plus[j] <- plus[j] + h
    minus[j] <- minus[j] - h
    H[, j] <- (score(plus) - score(minus)) / (2 * h)
  }
  info <- -0.5 * (H + t(H))
  ev <- eigen(info, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= max(ev) * 1e-8 || min(ev) <= 0)
    stop("Observed information is singular or nearly singular; standard errors are not reliable.")
  V <- solve(info)
  # Jacobian: from scaled (mu, tau2, logits) to original (mu, tau2, pi).
  J <- matrix(0, 3L * K, p)
  J[cbind(seq_len(K), seq_len(K))] <- scale
  J[cbind(K + seq_len(K), K + seq_len(K))] <- scale^2
  if (K > 1L) {
    for (k in seq_len(K)) for (j in seq_len(K - 1L))
      J[2L * K + k, 2L * K + j] <- pi[k] * (as.numeric(k == j) - pi[j])
  }
  Vout <- J %*% V %*% t(J)
  nm <- c(paste0("mu", seq_len(K)), paste0("tau2", seq_len(K)), paste0("pi", seq_len(K)))
  dimnames(Vout) <- list(nm, nm)
  se <- sqrt(pmax(0, diag(Vout)))
  names(se) <- nm
  estimates <- c(fit$mu, fit$tau2, fit$pi)
  names(estimates) <- nm
  list(se_mu = unname(se[seq_len(K)]),
       se_tau2 = unname(se[K + seq_len(K)]),
       se_pi = unname(se[2L * K + seq_len(K)]),
       estimates = estimates, se = se,
       vcov = Vout, information = info,
       min_information_eigenvalue = min(ev),
       method = "Inverse observed information; central differences of analytic score")
}


#' @export
table_fit.remix_cpp <- function(x, ...) {

  n <- x$n
  K <- x$K
  p <- x$npar
  LL <- x$logLik
  class <- x$classification
  # Posterior class probabilities
  posterior <- as.matrix(x$posterior)

  # Classification entropy
  # Normalized to [0, 1], with 1 indicating perfect classification
  if (K == 1L) {
    entropy <- 1
    #avg.mostlikely <- matrix(1L, dimnames = list("assigned.1", "meanprob.1"))
  } else {
    terms <- posterior * log(posterior)
    terms[!is.finite(terms)] <- 0
    entropy <- 1 + sum(terms) / (n * log(K))
  }

  sum.posterior <- colSums(posterior)
  sum.posterior <- data.frame(class = paste0("class", seq_len(K)), count = sum.posterior, proportion = sum.posterior/sum(sum.posterior))

  sum.mostlikely <- tabulate(x$classification, nbins = K)
  sum.mostlikely <- data.frame(class = paste0("class", seq_len(K)), count = sum.mostlikely, proportion = sum.mostlikely/sum(sum.mostlikely))

  individual <- cbind(posterior, class)
  colnames(individual) <- c(paste0("class", seq_len(K)), "predicted")

  # Classification certainty
  prob_max_case <- apply(posterior, 1, max)


  #avg.mostlikely <- avgprobs_mostlikely(posterior, class)
  avg.mostlikely <- t(sapply(1:ncol(posterior), function(i) {
    colMeans(posterior[class == i, , drop = FALSE])
  }))
  rownames(avg.mostlikely) <- paste0("assigned.", seq_len(nrow(avg.mostlikely)))
  colnames(avg.mostlikely) <- paste0("meanprob.", seq_len(ncol(avg.mostlikely)))
  avg.mostlikely[is.na(avg.mostlikely)] <- 0

  class_counts <- sum.mostlikely$count
  mostlikely.class <- avg.mostlikely * class_counts
  mostlikely.class <- mostlikely.class/matrix(colSums(avg.mostlikely * class_counts), ncol = ncol(mostlikely.class),
                      nrow = nrow(mostlikely.class), byrow = TRUE)
  rownames(mostlikely.class) <- paste0("assigned.", 1:nrow(mostlikely.class))
  colnames(mostlikely.class) <- paste0("avgprob.", 1:nrow(mostlikely.class))
  mostlikely.class <- t(mostlikely.class)


  # Sample-size adjusted BIC
  saBIC <- -2 * LL + p * log((n + 2) / 24)

  out <- data.frame(
    Minus2LogLikelihood = -2 * LL,
    n = n,
    Parameters = p,
    observedStatistics = 0,
    df = -p,
    RMSEASquared = 0,
    RMSEANull = 0.05,
    modelName = paste0("RE_Mixture", K),
    AIC = x$AIC,
    BIC = x$BIC,
    saBIC = saBIC,
    Classes = K,
    Entropy = entropy,
    prob_min = min(diag(mostlikely.class)),
    prob_max = max(diag(mostlikely.class)),
    n_min = min(sum.posterior$proportion),
    n_max = max(sum.posterior$proportion),
    LL = LL
  )
  out
}
