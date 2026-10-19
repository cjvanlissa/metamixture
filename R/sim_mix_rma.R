#' Simulate Data from a Random-Effects Mixture Model
#'
#' Simulates effect sizes from a finite mixture of normal random-effects
#' models, with class-specific means and residual heterogeneity.
#' @param n Number of observations to simulate.
#' @param vi_range Numeric vector specifying the minimum and maximum
#' sampling variances.
#' @param classes Number of mixture classes.
#' @param p_ratio Balance of class proportions. Default `1L` equally balances
#' class proportions; higher values assign a higher proportion to
#' the first class; e.g., `2L` means the first class will have twice the
#' proportion of cases as other classes.
#' @param r2 Explained variance of the true class difference (excluding sampling
#' variance `vi`). If `r2` is specified, `means` and `tau2` are ignored.
#' @param means Vector of class-specific mean effect sizes.
#' @param tau2 Vector of class-specific heterogeneity variances.
#' @return A data frame containing the simulated effect sizes (`yi`),
#'   sampling variances (`vi`), and true class memberships (`C`).
#' @examples
#' set.seed(123)
#' dat <- sim_mix_rma(n = 200, classes = 2,
#'                    p_ratio = 2,
#'                    means = c(0.2, 0.8),
#'                    tau2 = c(0.05, 0.10))
#' head(dat)
#' dat <- sim_mix_rma(n = 200, classes = 2,
#'                    p_ratio = 2,
#'                    means = c(0, 3),
#'                    tau2 = c(0.5, 2))
#' @importFrom stats var
#' @export
sim_mix_rma <- function(n = 100L, vi_range = c(.01, .10), classes = 2L, p_ratio = 1L, r2 = NULL, means = seq(from = 0, by = .2, length.out = classes), tau2 = rep(.1, classes)){
  ps <- c(p_ratio, rep(1, classes-1L))
  ps <- ps/sum(ps)
  C <- sample(seq_len(classes), size = n, replace = TRUE, prob = ps)
  if(!is.null(r2)){
    means <- 0:(classes-1L)
  }
  f_x <- (means[C])
  if(is.null(r2)){
    sigma_error <- sqrt(tau2)[C]
  } else {
    sigma_error <- sqrt(stats::var(f_x) * (1 - r2) / r2)
  }
  vi <- runif(n, min = vi_range[1], max = vi_range[2])
  y <- f_x + rnorm(n, sd = (sigma_error + sqrt(vi)))
  data.frame(yi = y, vi = vi, C = C)
}


# sim_mix_rma <- function(n = 100, vi_range = c(.01, .10), classes = 3, probs = rep(1/classes, classes), means = seq(from = 0, by = .2, length.out = classes), tau2 = rep(.1, classes)){
#   vi <- runif(n, vi_range[1], vi_range[2])
#   C <- sample.int(classes, n, replace = TRUE, prob = probs)
#   yi <- rnorm(n,
#               mean = means[C],
#               sd = sqrt(vi + tau2[C]))
#   data.frame(yi = yi, vi = vi, C = C)
# }
