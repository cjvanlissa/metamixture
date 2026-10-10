#' Simulate Data from a Random-Effects Mixture Model
#'
#' Simulates effect sizes from a finite mixture of normal random-effects
#' models, with class-specific means and residual heterogeneity.
#' @param n Number of observations to simulate.
#' @param vi_range Numeric vector specifying the minimum and maximum
#' sampling variances.
#' @param classes Number of mixture classes.
#' @param probs Vector of class membership probabilities.
#' @param means Vector of class-specific mean effect sizes.
#' @param tau2 Vector of class-specific heterogeneity variances.
#' @return A data frame containing the simulated effect sizes (`yi`),
#'   sampling variances (`vi`), and true class memberships (`C`).
#' @examples
#' set.seed(123)
#' dat <- sim_mix_rma(n = 200, classes = 2,
#'                    probs = c(0.7, 0.3),
#'                    means = c(0.2, 0.8),
#'                    tau2 = c(0.05, 0.10))
#' head(dat)
#'
#' @export
sim_mix_rma <- function(n = 100, vi_range = c(.01, .10), classes = 3, probs = rep(1/classes, classes), means = seq(from = 0, by = .2, length.out = classes), tau2 = rep(.1, classes)){
  vi <- runif(n, vi_range[1], vi_range[2])
  C <- sample.int(classes, n, replace = TRUE, prob = probs)
  yi <- rnorm(n,
              mean = means[C],
              sd = sqrt(vi + tau2[C]))
  data.frame(yi = yi, vi = vi, C = C)
}
