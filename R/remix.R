#' Random-Effects Mixture Meta-Analysis
#'
#' Fits finite mixture models for random-effects meta-analysis
#' using maximum likelihood estimation in `OpenMx`. Each mixture class
#' has its own mean effect size and residual heterogeneity variance.
#' Note that `OpenMx` must be loaded using `library(OpenMx)` before
#' calling this function.
#'
#' The estimation procedure is described in
#' Van Lissa, Garnier-Villarreal, & Anadria (2023): starting values
#' are determined by way of K-means clustering, then the model is
#' estimated using simulated annealing, followed by maximum
#' likelihood optimization.
#'
#' @param yi Numeric vector of observed effect sizes.
#' @param vi Numeric vector of known sampling variances.
#' @param classes Integer vector specifying the number of mixture
#' classes to fit. Defaults to \code{1L:2L}.
#' @param mu Whether component means are `varying` or `equal`.
#' @param tau2 Whether component heterogeneity variances are `varying` or
#' `equal`.
#' @param run Logical. If \code{TRUE} (the default), the model is
#' estimated using OpenMx. If \code{FALSE}, the unestimated
#' OpenMx model object is returned, allowing further modification
#' before estimation.
#' @return A fitted model object, or a `"mixture_list"` when multiple class
#' counts are specified.
#' @examples
#' \dontrun{
#' library(OpenMx)
#'
#' set.seed(7894)
#' df <- sim_mix_rma2(n = 200, classes = 2, p_ratio = 1, r2 = 1)
#'
#' # Conventional random-effects meta-analysis
#' res1 <- remix(df$yi, df$vi, classes = 1)
#' if(requireNamespace("metafor")){
#'   res_mf <- metafor::rma(yi = df$yi, vi = df$vi)
#'   res_mf$b[1,1] - res1$M$values[1, 1]
#'   res_mf$tau2 - res1$S$values[1, 1]
#' }
#'
#' # Two-class random-effects mixture model
#' res2 <- remix(df$yi, df$vi, classes = 2)
#' table_fit(res2)
#' table_results(res2)
#tmp <- metamixture:::random_starts(res2)
#table_fit(tmp)
#table_results(tmp)
#' }
#'
#' @references Van Lissa, C. J., Garnier-Villarreal, M., & Anadria, D. (2023).
#' Recommended Practices in Latent Class Analysis using the Open-Source
#' R-Package tidySEM. Structural Equation Modeling.
#' \doi{10.1080/10705511.2023.2250920}
#' @export
remix <- function(yi, vi,
                  classes = 2L,
                  mu = c("varying", "equal"),
                  tau2 = c("varying", "equal"),
                  run = TRUE) {
  intervals = TRUE
  mu <- match.arg(mu)
  tau2 <- match.arg(tau2)

  if (!requireNamespace("OpenMx", quietly = TRUE)) {
    stop("Package 'OpenMx' is required.")
  }
  s <- search()
  if(!any(grepl("OpenMx", s, fixed = TRUE))){
    message("Please run `library(OpenMx) before running `remix()`.")
    return(NULL)
  }
  if(length(classes) > 1){
    out <- lapply(classes, function(i){
      remix_mx(yi = yi, vi = vi, classes = i, mu = mu, tau2 = tau2, run = run)
    })
    attr(out, "tidySEM") <- "list"
    class(out) <- c("mixture_list", class(out))
    return(out)
  } else {
    # Input checks
    yi <- as.numeric(yi)
    vi <- as.numeric(vi)

    if (length(yi) != length(vi) ||
        length(yi) < 2L ||
        any(!is.finite(yi)) ||
        any(!is.finite(vi)) ||
        any(vi <= 0)) {
      stop("yi and vi must have equal lengths, non-missing, finite values, and vi > 0.")
    }

    if (length(classes) != 1L ||
        !is.finite(classes) ||
        classes < 1 ||
        classes != as.integer(classes)) {
      stop("'classes' must be a positive integer.")
    }

    classes <- as.integer(classes)
    n <- length(yi)

    if (classes > 1L){
      if(mu == "equal" && tau2 == "equal") {
        stop("For multiple classes, mu and tau2 cannot both be equal as this renders the mixture probabilities unidentified.")
      }
      if(length(unique(yi)) == 1L){
        stop("All effect sizes `yi` are identical. Cannot specify a mixture model.")
      }
    }



    # Observed data and known sampling variances
    objects <- list(
      OpenMx::mxMatrix(
        type = "Full", nrow = n, ncol = 1,
        free = FALSE, values = yi, name = "Y"
      ),
      OpenMx::mxMatrix(
        type = "Full", nrow = n, ncol = 1,
        free = FALSE, values = vi, name = "V"
      )
    )

    max_tau2 <- max(c(1e-7, diff(range(yi))^2))
    # Class-specific parameters and log-density algebras
    objects <- c(objects, list(
      OpenMx::mxMatrix(
        type = "Full",
        nrow = 1,
        ncol = classes,
        free = TRUE,
        values = 0,
        lbound = min(yi),
        ubound = max(yi),
        labels = if (mu == "varying") {
          paste0("mu", seq_len(classes))
        } else {
          rep("mu", classes)
        },
        name = "M"
      ),
      OpenMx::mxMatrix(
        type = "Diag",
        nrow = classes,
        ncol = classes,
        free = TRUE,
        values = min(max(1e-7, stats::var(yi) - mean(vi)), max_tau2),
        lbound = 1e-7,
        ubound = max_tau2,
        labels = if (tau2 == "varying") {
          paste0("tau2_", seq_len(classes))
        } else {
          rep("tau2", classes)
        },
        name = "S"
      )
    ))

    for (k in seq_len(classes)) {

      var_name <- paste0("Var", k)
      logd_name <- paste0("LogD", k)

      objects <- c(objects, list(

        OpenMx::mxAlgebraFromString(
          paste0("V + S[", k, ",", k, "]"),
          name = var_name
        ),

        OpenMx::mxAlgebraFromString(
          paste0(
            "-0.5 * (log(2 * pi * ", var_name, ") + ",
            "(Y - M[1,", k, "])^2 / ", var_name, ")"
          ),
          name = logd_name
        )

      ))
    }

    # Mixing proportions and likelihood
    if (classes == 1L) {

      # No free mixing proportion
      objects <- c(objects, list(
        OpenMx::mxMatrix(
          type = "Full", nrow = 1, ncol = 1,
          free = FALSE, values = 1,
          name = "weights"
        ),
        OpenMx::mxAlgebraFromString(
          "-2 * sum(LogD1)",
          name = "Minus2LL"
        )
      ))

    } else {

      # classes-1 free logit parameters
      objects <- c(objects, list(
        OpenMx::mxMatrix(
          type = "Full",
          nrow = 1, ncol = classes - 1L,
          free = TRUE,
          values = 0,
          lbound = -10,
          ubound = 10,
          labels = paste0("eta", seq_len(classes - 1L)),
          name = "Eta"
        )
      ))

      # Stable softmax: subtract maximum logit
      # The reference logit is zero.
      eta_expr <- paste0(
        "cbind(Eta, 0)"
      )

      objects <- c(objects, list(
        OpenMx::mxAlgebraFromString(
          eta_expr, name = "Logits"
        ),
        OpenMx::mxAlgebraFromString(
          "exp(Logits - max(Logits)) / sum(exp(Logits - max(Logits)))",
          name = "weights"
        )
      ))


      # Log of weighted component densities
      terms <- vapply(seq_len(classes), function(k) {
        paste0(
          "LogD", k, " + log(weights[1,", k, "])"
        )
      }, character(1))

      objects <- c(objects, list(
        OpenMx::mxAlgebraFromString(
          paste0("cbind(", paste(terms, collapse = ", "), ")"),
          name = "LogTerms"
        )
      ))

      # Compute row-wise maximum in separate algebra steps
      previous <- terms[1]

      for (k in 2:classes) {
        step_name <- paste0("MaxStep", k)

        expr <- paste0(
          "((", previous, ") + (", terms[k], ") + ",
          "abs((", previous, ") - (", terms[k], ")))/2"
        )

        objects <- c(objects, list(
          OpenMx::mxAlgebraFromString(
            expr,
            name = step_name
          )
        ))

        previous <- step_name
      }

      objects <- c(objects, list(
        OpenMx::mxAlgebraFromString(
          previous,
          name = "RowMax"
        ),

        # Repeat row-wise maximum across classes
        OpenMx::mxMatrix(
          type = "Unit",
          nrow = 1,
          ncol = classes,
          name = "Ones"
        ),

        OpenMx::mxAlgebra(
          RowMax %*% Ones,
          name = "RowMaxMatrix"
        ),

        # Exponentiate centered log densities
        OpenMx::mxAlgebra(
          exp(LogTerms - RowMaxMatrix),
          name = "ExpTerms"
        ),

        # Sum exponentials across classes
        OpenMx::mxMatrix(
          type = "Unit",
          nrow = classes,
          ncol = 1,
          name = "OnesCol"
        ),

        OpenMx::mxAlgebra(
          ExpTerms %*% OnesCol,
          name = "RowSumExp"
        ),

        # Final -2 log likelihood
        OpenMx::mxAlgebra(
          -2 * sum(RowMax + log(RowSumExp)),
          name = "Minus2LL"
        )
      ))
    }

    # Fit function
    objects <- c(objects, list(
      OpenMx::mxFitFunctionAlgebra("Minus2LL")
    ))

    # Confidence intervals
    if (intervals) {
      ci_names <- c("M", "S")
      if (classes > 1L) {
        ci_names <- c(ci_names, "weights")
      }
      objects <- c(objects, list(
        OpenMx::mxCI(ci_names)
      ))
    }

    # Construct model
    model <- do.call(
      OpenMx::mxModel,
      c(list(paste0("RE_Mixture", classes)), objects, OpenMx::mxData(data.frame(yi), type = "raw"))
    )
    # Improve starting values
    if(classes > 1){
      strts <- try(remix_starts(model), silent = TRUE)
      if(!inherits(strts, what = "try-error")){
        model <- strts
      } else {
        message("Could not improve starting values.")
      }
    }

    if (run) {
      model <- mx_run_remix(model)
    }
    attr(model, "tidySEM") <- "mixture"
    return(model)
  }
}
