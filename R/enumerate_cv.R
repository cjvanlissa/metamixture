remix_nlpd <- function(model, yi, vi) {

  logdens <- class_densities(model, yi, vi, normalize = FALSE)
  rowmax <- apply(logdens, 1, max)
  -sum(rowmax + log(rowSums(exp(logdens - rowmax))))
}


#' Cross-Validated Class Enumeration
#'
#' Selects the number of latent classes in random-effects mixture
#' meta-analysis using K-fold cross-validation. Candidate models
#' are compared based on their out-of-sample negative log predictive
#' density (NLPD).
#'
#' @param x An object for which a method exists, typically a `mixture_list` of
#' at least two fitted `remix` models with different numbers of latent classes.
#' The models must be based on the same observations and use the same
#' constraints on class-specific means and heterogeneity variances.
#' Models should be ordered by increasing number of classes,
#' starting with a one-class model.
#' @param folds Integer. Number of cross-validation folds.
#' Defaults to 5. Must be between 2 and the number of observations.
#' @details
#' The function randomly partitions the observations into `folds`
#' approximately equal-sized subsets. For each candidate number of
#' classes, the model is re-estimated on the training observations
#' and evaluated on the held-out observations using the negative log
#' predictive density. Smaller values indicate better out-of-sample predictive
#' performance.
#'
#' The negative log predictive densities are summed across folds
#' for each candidate model. The model with the smallest total
#' score is selected.
#'
#' The mean and heterogeneity variance constraints are inferred
#' from the second model in `x` and applied consistently to all
#' models fitted during cross-validation.
#'
#' Because fold assignments are random, results may differ between
#' calls. Use [set.seed()] before calling the function to obtain
#' reproducible fold assignments and random starting values.
#'
#' @return
#' A list containing:
#' \describe{
#' \item{best_class}{Integer index of the candidate model with
#'   the lowest total cross-validated negative log predictive
#'   density. If models are ordered from one to the maximum
#'   number of classes, this index corresponds to the selected
#'   number of classes.}
#' \item{scores}{Named numeric vector containing the total
#'   negative log predictive density across folds for each
#'   candidate model.}
#' \item{fold_scores}{Numeric matrix containing the negative
#'   log predictive density for each fold (rows) and candidate
#'   model (columns).}
#' \item{fold_id}{Integer vector indicating the fold assigned
#'   to each observation.}
#' }
#'
#' @examples
#' \dontrun{
#' library(metamixture)
#'
#' set.seed(123)
#'
#' # Simulate meta-analytic data
#' dat <- sim_remix_fixed_tau(n = 100, es = 2)
#'
#' # Estimate candidate models
#' models <- lapply(1:3, function(k) {
#' remix(
#'   yi = dat$yi,
#'   vi = dat$vi,
#'   classes = k,
#'   tau2 = "equal"
#' )
#' })
#'
#' # Five-fold cross-validation
#' cv <- cv_enumeration(models, folds = 5)
#'
#' # Selected number of classes
#' cv$best_class
#'
#' # Total predictive scores
#' cv$scores
#'
#' # Scores by fold and candidate model
#' cv$fold_scores
#' }
#'
#' @export
cv_enumeration <- function(x, folds = 5L) {
  UseMethod("cv_enumeration", x)
}

#' @method cv_enumeration mixture_list
#' @export
cv_enumeration.mixture_list <- function(x, folds = 5L) {

  if (length(x) < 2L) {
    stop("At least two candidate models are required.")
  }

  yi <- as.numeric(x[[1]]$Y$values)
  vi <- as.numeric(x[[1]]$V$values)
  n <- length(yi)
  if (folds < 2L || folds > n) {
    stop("'folds' must be between 2 and n.")
  }
  mu_equal <- length(unique(x[[2]]$M$labels[1,])) == 1L
  tau_equal <- length(unique(diag(x[[2]]$S$labels))) == 1L

  fold_id <- sample(rep(seq_len(folds), length.out = n))

  cv_scores <- vapply(seq_along(x), function(K) {

    vapply(seq_len(folds), function(f) {

      train <- fold_id != f
      test <- !train

      fit <- remix(yi = yi[train], vi = vi[train], classes = K, mu = c("varying", "equal")[mu_equal+1L], tau2 = c("varying", "equal")[tau_equal+1L], random_starts = FALSE)
      if(!mx_diagnostics(fit, check = "fit")$checks_pass){
        fit <- tidySEM::random_starts(fit, verbose = FALSE)
      }

      # Evaluate held-out predictive score
      remix_nlpd(
        fit,
        yi = yi[test],
        vi = vi[test]
      )

    }, numeric(1))

  }, numeric(folds))

  # Total out-of-sample negative log predictive density
  total_scores <- colSums(cv_scores)

  names(total_scores) <- tryCatch(sapply(x, function(m) m@name), error = function(e) paste0("model ", seq_along(x)))

  list(
    best_class = which.min(total_scores),
    scores = total_scores,
    fold_scores = cv_scores,
    fold_id = fold_id
  )
}
