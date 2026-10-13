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
#' @param method Estimation method: `"cpp"` or `"OpenMx"`.
#' @param run Logical. If \code{TRUE} (the default), the model is
#' estimated using OpenMx. If \code{FALSE}, the unestimated
#' OpenMx model object is returned, allowing further modification
#' before estimation. Ignored for C++.
#' @param options Named list of arguments passed to `remix_cpp()`.
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
#'   res_mf$b[1,1] - res1$Mu1$values[1,1]
#'   res_mf$tau2 - res1$Tau1$values[1,1]
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
remix <- function(yi, vi, classes = 2L,
                  method = c("cpp", "OpenMx"),
                  run = TRUE,
                  options = list()) {

  method <- match.arg(method)
  classes <- as.integer(classes)

  if (method == "OpenMx") {

    if (length(options)) {
      stop("'options' is only supported for `method = 'cpp'`.")
    }

    return(remix_mx(
      yi = yi,
      vi = vi,
      classes = classes,
      run = run
    ))
  }

  # C++ method
  valid_options <- c(
    "nstart", "maxit", "tol", "grid_size",
    "profile_tol", "joint", "joint_maxit",
    "joint_factr", "joint_pgtol", "init",
    "keep_trace", "compute_se", "se_rel_step",
    "se_boundary_tol"
  )

  unknown <- setdiff(names(options), valid_options)

  if (length(unknown)) {
    stop("Unknown C++ options: ",
         paste(unknown, collapse = ", "))
  }

  fit_one <- function(k) {
    do.call(
      remix_cpp,
      c(list(yi = yi, vi = vi, K = k), options)
    )
  }

  if (length(classes) == 1L) {
    return(fit_one(classes))
  }

  out <- lapply(classes, fit_one)

  names(out) <- paste0("RE_Mixture", classes)

  attr(out, "tidySEM") <- "list"
  class(out) <- c("mixture_list", "list")
  out
}
