#' Cluster-Robust Standard Errors for `remix` Models
#'
#' Computes cluster-robust sandwich standard errors and confidence intervals
#' for random-effects mixture models fitted using [remix()].
#' This allows for dependence between effect sizes belonging to the same
#' cluster, such as multiple effect sizes from the same study.
#
#' @param model A fitted OpenMx model created using [remix()].
#' @param cluster Vector identifying the cluster membership of each effect
#' size, e.g., the study the effect size came from. Must have the same length
#' and order as the observed effect sizes.
#' @param type Cluster-robust covariance estimator. Either `"CR1"`
#' (default), applying the correction G/(G-1), or `"CR0"`, the
#' unadjusted sandwich estimator, where G is the number of clusters (or,
#' studies).
#' @param rel_step Relative step size for numerical differentiation of
#' the analytic score function. Defaults to `1e-4`.
#' @return An `MxModel` with `attr("se_robust")` containing the standard errors,
#' and `attr("sandwich")` with
#' `vcov`, the cluster-robust covariance matrix, `vcov_free`, the cluster-robust
#' covariance matrix for the unique free model parameters, `information`,
#' observed Fisher information matrix (negative Hessian of the log-likelihood),
#' `cluster_scores`, the matrix of first derivatives of the log-likelihood for
#' each cluster (rows) by free parameters (columns); `n` number of observed
#' effect sizes, `clusters` number clusters G, `npar`, number of freely
#' estimated parameters, `type` passed from input.
#' @details
#' Note that cluster-robust standard errors provide consistent inference for
#' dependent effect sizes but do not affect the fitted parameter estimates, nor
#' model the dependency in the data.
#' CR0 and CR1 are large-sample approximations and may be unreliable when the
#' number of clusters is small or mixture parameters are weakly identified.
#'
#' @examples
#' \dontrun{
#' library(OpenMx)
#' library(metamixture)
#'
#' set.seed(1)
#'
#' n <- 100
#' df <- metamixture::sim_remix_fixed_tau(n, es = 2)
#'
#' res1 <- remix(
#'   yi = df$yi,
#'   vi = df$vi,
#'   classes = 2,
#'   tau2 = "equal"
#' )
#'
#' # Original sample: one observation per cluster
#' sw1 <- sandwich_remix(
#'   res1,
#'   cluster = seq_len(n),
#'   type = "CR0"
#' )
#'
#' # Duplicate observations exactly
#' df2 <- rbind(df, df)
#' df2$study <- rep(seq_len(n), 2)
#'
#' res2 <- remix(
#'   yi = df2$yi,
#'   vi = df2$vi,
#'   classes = 2,
#'   tau2 = "equal"
#' )
#'
#' sw2 <- sandwich_remix(
#'   res2,
#'   cluster = df2$study,
#'   type = "CR0"
#' )
#'
#' # Compare robust covariance matrices
#' all.equal(
#'   attr(sw1, "sandwich")$vcov,
#'   attr(sw2, "sandwich")$vcov,
#'   tolerance = 1e-3
#' )
#' }
#'
#' @export
sandwich_remix <- function(model, cluster,
                           type = c("CR1", "CR0"),
                           rel_step = 1e-4) {

  type <- match.arg(type)

  # Observed data
  yi <- as.numeric(model$Y$values)
  vi <- as.numeric(model$V$values)
  if(!length(vi) == length(cluster)){
    stop("Argument `cluster` must be a vector of the same length as the number of observations used to estimate `model`.")
  }
  n <- length(yi)
  K <- ncol(model$M$values)

  cluster <- as.integer(factor(cluster))
  G <- length(unique(cluster))

  mu_equal <- length(unique(model$M$labels[1, ])) == 1L
  tau_equal <- length(unique(diag(model$S$labels))) == 1L

  # Get parameters
  mu <- as.numeric(model$M$values[1, ])
  tau2 <- diag(model$S$values)

  pi <- as.numeric(OpenMx::mxEvalByName("weights", model))

  if (any(!is.finite(c(mu, tau2, pi))) || any(tau2 < 0) || any(pi <= 0)) {
    stop("Invalid fitted parameters.")
  }

  # Reduced parameterization
  # Common parameters appear only once.
  m <- if (mu_equal) mu[1] else mu
  t <- if (tau_equal) tau2[1] else tau2

  eta <- if (K > 1L) {
    log(pi[seq_len(K - 1L)] / pi[K])
  } else {
    numeric()
  }

  theta <- c(m, t, eta)

  n_mu <- length(m)
  n_tau <- length(t)
  p <- length(theta)

  # Convert reduced parameters to full parameters
  unpack <- function(th) {

    mu_par <- th[seq_len(n_mu)]

    tau_par <- th[n_mu + seq_len(n_tau)]

    eta_par <- if (K > 1L) {
      th[n_mu + n_tau + seq_len(K - 1L)]
    } else {
      numeric()
    }

    mu_full <- if (mu_equal) rep(mu_par, K) else mu_par
    tau_full <- if (tau_equal) rep(tau_par, K) else tau_par

    logits <- c(eta_par, 0)
    z <- exp(logits - max(logits))
    probs <- z / sum(z)

    list(mu = mu_full, tau2 = tau_full, pi = probs)
  }

  # Analytic observation-level score matrix
  # Dimensions: n observations x p free parameters
  score_obs <- function(th) {

    par <- unpack(th)

    mu <- par$mu
    tau2 <- par$tau2
    pi <- par$pi

    s <- outer(vi, tau2, "+")

    if (any(s <= 0)) {
      stop("Nonpositive total variance.")
    }

    resid <- outer(yi, mu, "-")

    logdens <- sweep(-0.5 * (log(2 * base::pi) + log(s) + resid^2 / s), 2, log(pi), "+")

    # Stable posterior probabilities
    rowmax <- apply(logdens, 1, max)

    a <- exp(logdens - rowmax)
    r <- a / rowSums(a)

    # Mean scores
    score_mu <- r * resid / s

    if (mu_equal) {
      score_mu <- matrix(rowSums(score_mu), ncol = 1)
    }

    # Heterogeneity variance scores
    score_tau <- 0.5 * r * (
      resid^2 / s^2 - 1 / s
    )

    if (tau_equal) {
      score_tau <- matrix(rowSums(score_tau), ncol = 1)
    }

    # Mixing proportion logit scores
    score_eta <- if (K > 1L) {
      sweep(
        r[, seq_len(K - 1L), drop = FALSE],
        2,
        pi[seq_len(K - 1L)],
        "-"
      )
    } else {
      matrix(numeric(0), nrow = n, ncol = 0)
    }

    cbind(score_mu, score_tau, score_eta)
  }

  # Total analytic score
  total_score <- function(th) {
    colSums(score_obs(th))
  }

  # Numerical Hessian of log-likelihood
  H <- matrix(0, p, p)

  for (j in seq_len(p)) {

    h <- rel_step * max(1, abs(theta[j]))

    # Keep tau2 perturbations positive
    if (j > n_mu && j <= n_mu + n_tau) {
      h <- min(h, theta[j] * 0.25)
    }

    if (!is.finite(h) || h <= 0) {
      stop("Cannot differentiate at a boundary estimate.")
    }

    plus <- minus <- theta

    plus[j] <- plus[j] + h
    minus[j] <- minus[j] - h

    H[, j] <- (
      total_score(plus) -
        total_score(minus)
    ) / (2 * h)
  }

  # Observed information
  info <- -(H + t(H)) / 2

  ev <- eigen(
    info,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  if (min(ev) <= 0 ||
      min(ev) <= max(ev) * 1e-8) {
    stop("Observed information is singular or nearly singular.")
  }

  bread <- solve(info)

  # Aggregate observation scores by cluster
  scores <- score_obs(theta)

  cluster_scores <- rowsum(
    scores,
    group = cluster,
    reorder = TRUE
  )

  # Sandwich topping
  hummus <- crossprod(cluster_scores)

  # CR0
  V <- bread %*% hummus %*% bread

  # CR1 correction
  if (type == "CR1") {
    V <- V * G / (G - 1)
  }

  # Jacobian from reduced parameters to
  # class-specific mu, tau2, and pi
  J <- matrix(0, 3L * K, p)

  for (k in seq_len(K)) {

    # Means
    j <- if (mu_equal) 1L else k
    J[k, j] <- 1

    # Heterogeneity variances
    j <- if (tau_equal) 1L else k
    J[K + k, n_mu + j] <- 1
  }

  # Softmax Jacobian
  if (K > 1L) {
    for (k in seq_len(K)) {
      for (j in seq_len(K - 1L)) {

        J[2L * K + k, n_mu + n_tau + j] <-
          pi[k] * (as.numeric(k == j) - pi[j])
      }
    }
  }

  # Covariance matrix on original parameter scale
  Vout <- J %*% V %*% t(J)
  Vout <- (Vout + t(Vout)) / 2

  nm <- c(
    paste0("mu", seq_len(K)),
    paste0("tau2_", seq_len(K)),
    paste0("weights[1,", seq_len(K), "]")
  )

  dimnames(Vout) <- list(nm, nm)

  se <- sqrt(pmax(0, diag(Vout)))
  names(se) <- nm
  ses <- data.frame(
    label = nm,
    se = se,
    type = "sandwich"
  )
  attr_list <- list(
    vcov = Vout,
    vcov_free = V,
    information = info,
    cluster_scores = cluster_scores,
    n = n,
    clusters = G,
    npar = p,
    type = type
  )

  attr(model, "se_robust") <- ses
  attr(model, "sandwich") <- attr_list
  return(model)
}
