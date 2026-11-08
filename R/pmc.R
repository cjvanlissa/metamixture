#' Class Enumeration Using Predictive Surprisal
#'
#' Compares random-effects mixture meta-analysis models using
#' simulation-based predictive surprisal. For each candidate model,
#' repeated datasets are simulated, and kernel density estimation
#' (KDE) is used to evaluate the surprisal of the observed effect
#' sizes under the simulated distribution.
#'
#' @param x A list of at least two fitted `remix` models, typically
#' representing different numbers of latent classes. All models
#' must be fitted to the same observed effect sizes.
#' @param ... Additional arguments passed to [surprisal()].
#' @param reps Integer. Number of simulated datasets generated
#' for each candidate model. Defaults to 100.
#' @param ci Numeric. Confidence level for the bootstrap confidence
#' intervals of differences in mean surprisal. Defaults to 0.95.
#' @param nboot Integer. Number of bootstrap replications used
#' to estimate confidence intervals for differences in mean
#' surprisal. Defaults to 100.
#'
#' @details
#' For each candidate model, the function generates `reps`
#' simulated datasets.
#'
#' For each simulated dataset, [surprisal()] estimates the
#' distribution of the simulated effect sizes using kernel
#' density estimation and evaluates the surprisal of the
#' observed effect sizes under this distribution.
#'
#' The mean surprisal across simulated datasets is used as
#' a measure of predictive discrepancy. Smaller values indicate
#' that the simulated distributions assign greater density to
#' the observed effect sizes and therefore provide a better
#' approximation to the observed data.
#'
#' Two types of model comparisons are performed:
#' \describe{
#' \item{dif_seq}{Compares each model with the preceding model
#'   in the supplied list.}
#' \item{dif_one}{Compares each model with the first model
#'   in the supplied list.}
#' }
#'
#' For each comparison, the difference is calculated as
#' the mean surprisal of the alternative model minus the
#' mean surprisal of the reference model:
#'
#' \deqn{\Delta = \bar{S}_{alt} - \bar{S}_{null}}
#'
#' Negative differences favor the alternative model.
#'
#' Confidence intervals are obtained using a nonparametric
#' percentile bootstrap. The simulated surprisal scores
#' for the two models are independently resampled with
#' replacement, and the difference in their means is
#' recalculated for each bootstrap replication.
#'
#' A comparison is marked with an asterisk (`*`) when the
#' upper confidence limit is below zero, indicating that
#' the alternative model has significantly lower mean
#' surprisal at the specified confidence level.
#'
#' If the `future.apply` package is installed, simulations
#' are evaluated using [future.apply::future_replicate()],
#' allowing parallel execution when an appropriate
#' `future` plan has been configured. Otherwise, the
#' function uses sequential replication.
#'
#' The function uses random simulation and bootstrap
#' resampling. Use [set.seed()] for reproducibility.
#'
#' @return
#' A data.frame of class `pmc_df` containing:
#' \describe{
#' \item{comparison}{Comparison type: `dif_seq` or `dif_one`.}
#' \item{null}{Name of the reference model.}
#' \item{alt}{Name of the alternative model.}
#' \item{null_stat}{Mean simulated-data surprisal for the
#'   reference model.}
#' \item{alt_stat}{Mean simulated-data surprisal for the
#'   alternative model.}
#' \item{difference}{Difference in mean surprisal
#'   (`alt_stat - null_stat`).}
#' \item{lb}{Lower bootstrap confidence limit.}
#' \item{ub}{Upper bootstrap confidence limit.}
#' \item{sig}{An asterisk (`*`) if the upper confidence
#'   limit is below zero; otherwise an empty string.}
#' }
#'
#' @examples
#' \dontrun{
#' library(metamixture)
#'
#' set.seed(123)
#'
#' # Simulate meta-analytic data
#' dat <- sim_remix_fixed_tau(k = 100, es = 2)
#'
#' # Fit candidate models
#' models <- lapply(1:3, function(k) {
#' remix(
#'   yi = dat$yi,
#'   vi = dat$vi,
#'   classes = k,
#'   tau2 = "equal"
#' )
#' })
#'
#' # Compare models using predictive surprisal
#' pmc <- pmc_surprisal(
#' models,
#' reps = 100,
#' nboot = 500
#' )
#'
#' # Model comparisons
#' pmc
#'
#' # Simulated surprisal scores
#' head(attr(pmc, "scores"))
#' }
#'
#'
#' @importFrom stats median quantile
#' @export
pmc_surprisal <- function(x, ..., reps = 100L, ci = 0.95,
                   nboot = 100L) {

  if (length(x) < 2L) {
    stop("At least two models are required.")
  }

  if (length(reps) != 1L || !is.finite(reps) ||
      reps < 2L || reps != as.integer(reps)) {
    stop("'reps' must be an integer >= 2.")
  }

  if (length(ci) != 1L || !is.finite(ci) ||
      ci <= 0 || ci >= 1) {
    stop("'ci' must be between 0 and 1.")
  }

  if (length(nboot) != 1L || !is.finite(nboot) ||
      nboot < 2L || nboot != as.integer(nboot)) {
    stop("'nboot' must be an integer >= 2.")
  }

  K <- length(x)

  nams <- tryCatch(
    vapply(x, function(m) m@name, character(1)),
    error = function(e) paste0("model ", seq_len(K))
  )

  if (anyDuplicated(nams)) {
    nams <- make.unique(nams)
  }

  # Observed effect sizes
  yi <- as.numeric(x[[1]]$Y$values)

  # One simulated score
  one_rep <- function(model) {

    sim <- simulate_remix(
      model,
      add_vi = TRUE
    )

    # Density estimated from simulated observations,
    # evaluated at observed effects
    score <- surprisal(
      x = sim$yi,
      y = yi,
      ...
    )

    if (length(score) != 1L ||
        !is.finite(score)) {
      stop("Surprisal must return one finite value.")
    }

    as.numeric(score)
  }

  # Simulate scores for each model
  rep_stat <- vapply(seq_along(x), function(i) {

    model <- x[[i]]

    if (requireNamespace("future.apply", quietly = TRUE)) {

      scores <- future.apply::future_replicate(
        n = reps,
        expr = one_rep(model),
        future.seed = TRUE,
        simplify = FALSE
      )

      unlist(scores, use.names = FALSE)

    } else {

      replicate(reps, one_rep(model))

    }

  }, FUN.VALUE = numeric(reps))

  colnames(rep_stat) <- nams

  # Mean expected surprisal
  means <- colMeans(rep_stat)
  names(means) <- nams

  # Comparisons:
  # successive models and each model against the first
  comparisons <- rbind(
    data.frame(
      comparison = "dif_seq",
      null_idx = seq_len(K - 1L),
      alt_idx = seq.int(2L, K)
    ),
    data.frame(
      comparison = "dif_one",
      null_idx = rep(1L, K - 1L),
      alt_idx = seq.int(2L, K)
    )
  )

  # Bootstrap CI for difference in mean scores
  compare_scores <- function(a, b) {

    delta <- mean(b) - mean(a)

    boot <- replicate(nboot, {

      aa <- sample(a, length(a), replace = TRUE)
      bb <- sample(b, length(b), replace = TRUE)

      mean(bb) - mean(aa)
    })

    limits <- stats::quantile(
      boot,
      probs = c((1 - ci) / 2, 1 - (1 - ci) / 2),
      names = FALSE
    )

    c(
      difference = delta,
      lb = limits[1],
      ub = limits[2]
    )
  }

  comparisons_stats <- t(vapply(
    seq_len(nrow(comparisons)),
    function(j) {

      compare_scores(
        rep_stat[, comparisons$null_idx[j]],
        rep_stat[, comparisons$alt_idx[j]]
      )

    },
    FUN.VALUE = c(
      difference = 0,
      lb = 0,
      ub = 0
    )
  ))

  out <- data.frame(
    comparison = comparisons$comparison,
    null = nams[comparisons$null_idx],
    alt = nams[comparisons$alt_idx],
    null_stat = means[comparisons$null_idx],
    alt_stat = means[comparisons$alt_idx],
    difference = comparisons_stats[, "difference"],
    lb = comparisons_stats[, "lb"],
    ub = comparisons_stats[, "ub"],
    sig = ifelse(
      comparisons_stats[, "ub"] < 0,
      "*",
      ""
    ),
    row.names = NULL
  )

  class(out) <- c("pmc_df", class(out))
  attr(out, "stat") <- "mean simulated-data KDE surprisal"
  attr(out, "scores") <- rep_stat
  out
}
