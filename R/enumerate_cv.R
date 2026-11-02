# remix_surprisal <- function(model, yi, vi) {
#
#   mu <- as.numeric(OpenMx::mxEval(M, model))
#   tau2 <- diag(OpenMx::mxEval(S, model))
#   pi <- as.numeric(OpenMx::mxEval(weights, model))
#
#   K <- length(pi)
#
#   logdens <- vapply(seq_len(K), function(k) {
#     log(pi[k]) +
#       dnorm(
#         yi,
#         mean = mu[k],
#         sd = sqrt(vi + tau2[k]),
#         log = TRUE
#       )
#   }, numeric(length(yi)))
#
#   logdens <- matrix(logdens, ncol = K)
#
#   rowmax <- apply(logdens, 1, max)
#
#   -sum(
#     rowmax +
#       log(rowSums(exp(logdens - rowmax)))
#   )
# }
#
#
# cv_enumeration <- function(x, folds = 5){
#
#   yi <- x[[1]]$Y$values
#   vi <- x[[1]]$V$values
#   fold_id <- sample.int(folds, size = length(yi), replace = TRUE)
#
#
# cv_scores <- sapply(seq_along(x), function(K){
#   sapply(seq_len(folds), function(f){
#
#     train <- fold_id != f
#     test <- !train
#     n_cv <- sum(train)
#     mat_Y <-  OpenMx::mxMatrix(
#       type = "Full", nrow = n_cv, ncol = 1,
#       free = FALSE, values = yi[train], name = "Y"
#     )
#     mat_V <- OpenMx::mxMatrix(
#       type = "Full", nrow = n_cv, ncol = 1,
#       free = FALSE, values = vi[train], name = "V"
#     )
#     mat_data <- OpenMx::mxData(data.frame(yi = yi[train]), type = "raw")
#     fit <- mxModel(x[[K]], mat_Y, mat_V, mat_data)
#     fit <- mxRun(fit, silent = TRUE, suppressWarnings = TRUE)
#     remix_surprisal(fit, yi[test], vi[test])
#   })
# })
#
# which.min((colSums(cv_scores)+ (log(x[[1]]@data$numObs) * sapply(seq_along(x), function(i)length(x[[i]]@output$estimate)))))
# }
#



remix_nlpd <- function(model, yi, vi) {

  mu <- as.numeric(OpenMx::mxEvalByName("M", model))
  tau2 <- diag(OpenMx::mxEvalByName("S", model))
  pi <- as.numeric(OpenMx::mxEvalByName("weights", model))

  K <- length(pi)

  logdens <- vapply(seq_len(K), function(k) {
    log(pi[k]) +
      stats::dnorm(
        yi,
        mean = mu[k],
        sd = sqrt(vi + tau2[k]),
        log = TRUE
      )
  }, numeric(length(yi)))

  logdens <- matrix(logdens, nrow = length(yi), ncol = K)

  rowmax <- apply(logdens, 1, max)

  -sum(
    rowmax +
      log(rowSums(exp(logdens - rowmax)))
  )
}


cv_enumeration <- function(x, folds = 5L) {
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
      if(!mx_convergence(fit, check = "fit")$checks_pass){
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
