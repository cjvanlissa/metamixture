#' Simulate Data from a Random-Effects Mixture Model
#'
#' Simulates effect sizes from a finite mixture of normal random-effects
#' models, with class-specific means and residual heterogeneity.
#' @param k Number of effect sizes to simulate.
#' @param classes Number of mixture classes.
#' @param vi_range Numeric vector specifying the minimum and maximum
#' sampling variances.
#' @param p_ratio Numeric, balance of class proportions. Default `1` equally
#' balances class proportions; higher values assign a higher proportion to
#' the first class; e.g., `2` means the first class will have twice the
#' proportion of cases as other classes.
# @param r2 Explained variance of the true class difference (excluding sampling
# variance `vi`). If `r2` is specified, `means` and `tau2` are ignored.
# @param means Vector of class-specific mean effect sizes.
#' @param es Numeric, effect size (distance between class centroids).
# @param tau2 Vector of class-specific heterogeneity variances.
#' @return A data frame containing the simulated effect sizes (`yi`),
#'   sampling variances (`vi`), and true class memberships (`C`).
#' @examples
#' set.seed(123)
#' dat <- sim_remix_fixed_tau(k = 200, classes = 2, p_ratio = 2, means = c(0.2, 0.8),
#'                    tau2 = c(0.05, 0.10))
#' head(dat)
#' @export
sim_remix_fixed_tau <- function(k = 100L, classes = 2L, p_ratio = 1, vi_range = c(.001, .01), es = 1) {
  if(classes > 1){
    ps <- c(p_ratio, rep(1, classes - 1L))
    ps <- ps / sum(ps)
  } else {
    ps <- 1
  }
  mu <- (seq_len(classes) - 1L) * es
  vi_range <- vi_range
  C <- sample(seq_len(classes), size = k, replace = TRUE, prob = ps)

  vi_indiv <- runif(k, min = vi_range[2], max = vi_range[2])
  yi <- rnorm(k, mean = mu[C], sd = sqrt(1 + vi_indiv))

  data.frame(
    yi = yi,
    vi = vi_indiv,
    C = C
  )
}


# sim_mix_rma <- function(n = 100, vi_range = c(.01, .10), classes = 3, probs = rep(1/classes, classes), means = seq(from = 0, by = .2, length.out = classes), tau2 = rep(.1, classes)){
#   vi <- runif(n, vi_range[1], vi_range[2])
#   C <- sample.int(classes, n, replace = TRUE, prob = probs)
#   yi <- rnorm(n,
#               mean = means[C],
#               sd = sqrt(vi + tau2[C]))
#   data.frame(yi = yi, vi = vi, C = C)
# }
