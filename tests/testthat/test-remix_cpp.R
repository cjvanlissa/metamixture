RNGkind(
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

simulate_re_mixture <- function(n, pi, mu, tau2, vi_range = c(0.005, 0.025)) {
  K <- length(pi)
  stopifnot(length(mu) == K, length(tau2) == K,
            all(pi > 0), abs(sum(pi) - 1) < 1e-10,
            all(tau2 >= 0))
  vi <- stats::runif(n, vi_range[1], vi_range[2])
  z <- sample.int(K, n, replace = TRUE, prob = pi)
  yi <- stats::rnorm(n, mean = mu[z], sd = sqrt(vi + tau2[z]))
  list(yi = yi, vi = vi, z = z)
}


test_that("one-class random-effects model recovers its parameters", {
  set.seed(828)
  dat <- simulate_re_mixture(100, pi = 1, mu = 0.35, tau2 = 0.12)
  fit <- metamixture:::remix_cpp(dat$yi, dat$vi, K = length(pi),
    init = list(mu = .35, tau2 = .12, pi = 1), compute_se = FALSE
  )
  expect_equal(fit$K, 1L)
  expect_equal(fit$npar, 2L)
  expect_equal(fit$pi, 1, tolerance = 1e-12)
  expect_true(abs(fit$mu-0.35) < 0.03)
  expect_true(abs(fit$tau2-.12) < 0.03)
  expect_equal(fit$posterior, matrix(1, nrow = length(dat$yi), ncol = 1))
  expect_true(is.finite(fit$logLik))
})

test_that("one-class model works without starts", {
  set.seed(73368)
  dat <- simulate_re_mixture(100, pi = 1, mu = 0.35, tau2 = 0.12)
  fit <- metamixture:::remix_cpp(dat$yi, dat$vi, K = length(pi), compute_se = FALSE)
  expect_equal(fit$K, 1L)
  expect_equal(fit$npar, 2L)
  expect_equal(fit$pi, 1, tolerance = 1e-12)
  expect_true(abs(fit$mu-0.35) < 0.04)
  expect_true(abs(fit$tau2-.12) < 0.03)
  expect_equal(fit$posterior, matrix(1, nrow = length(dat$yi), ncol = 1))
  expect_true(is.finite(fit$logLik))
})

test_that("two classes recover varying means and heterogeneity", {
  set.seed(831)
  dat <- simulate_re_mixture(n = 100, pi = c(0.60, 0.40), mu = c(-1.2, 1.4),
                             tau2 = c(0.035, 0.16))
  fit <- metamixture:::remix_cpp(dat$yi, dat$vi, K = 2, compute_se = FALSE)
  expect_equal(fit$npar, 5L)
  expect_equal(fit$mu, c(-1.2, 1.4), tolerance = 0.2)
  expect_true(all(abs(fit$tau2 - c(0.035, 0.16)) < .06))
  expect_true(all(abs(fit$pi- c(0.60, 0.40)) < c(.05, .05)))
  expect_equal(sum(fit$pi), 1, tolerance = 1e-10)
  expect_equal(rowSums(fit$posterior), rep(1, length(dat$yi)),
               tolerance = 1e-10)
  expect_equal(sum(fit$effective_n), length(dat$yi), tolerance = 1e-8)
  expect_equal(fit$AIC, -2 * fit$logLik + 2 * fit$npar,
               tolerance = 1e-8)
})

test_that("equal means and varying heterogeneity are recovered", {
  set.seed(14885)
  dat <- simulate_re_mixture(n = 100, pi = c(0.60, 0.40), mu = c(0.25, 0.25),
                             tau2 = c(0.015, 1.0), vi_range = c(0.002, 0.015))
  fit <- metamixture:::remix_cpp(dat$yi, dat$vi, K = 2, compute_se = FALSE, mu = "equal")
  expect_equal(fit$npar, 4L)
  expect_equal(fit$mu[1], fit$mu[2], tolerance = 1e-10)
  expect_true(abs(fit$mu[1] - .25) < 0.05)
  expect_true(all(abs(fit$tau2 - c(0.015, 1.0)) < 0.4))
  expect_equal(fit$pi, c(.6, .4), tolerance = 0.10)
})

test_that("varying means and equal heterogeneity are recovered", {
  set.seed(953)
  dat <- simulate_re_mixture(n = 100, pi = c(0.55, 0.45), mu = c(-1.0, 1.2),
  tau2 = c(0.09, 0.09))
  fit <- metamixture:::remix_cpp(dat$yi, dat$vi, K = 2, compute_se = FALSE, tau2 = "equal")
  expect_equal(fit$npar, 4L)
  expect_equal(fit$tau2[1], fit$tau2[2], tolerance = 1e-10)
  expect_equal(fit$mu, c(-1,1.2), tolerance = 0.1)
  expect_true(all(abs(fit$tau2 - c(0.09, 0.09)) < 0.02))
  expect_true(all(abs(fit$pi - c(0.55, 0.45)) < .06))
})

test_that("random starts respect R's random seed", {
  set.seed(67)
  dat <- simulate_re_mixture(35, c(0.6, 0.4), c(-1.0, 1.0),
                             c(0.05, 0.10))
  set.seed(777)
  a <- metamixture:::remix_cpp(dat$yi, dat$vi, K = 2, compute_se = FALSE)
  set.seed(777)
  b <- metamixture:::remix_cpp(dat$yi, dat$vi, K = 2, compute_se = FALSE)
  expect_equal(a$mu, b$mu, tolerance = 1e-12)
  expect_equal(a$tau2, b$tau2, tolerance = 1e-12)
  expect_equal(a$pi, b$pi, tolerance = 1e-12)
  expect_equal(a$starts$logLik, b$starts$logLik, tolerance = 1e-12)
})

test_that("unidentified and invalid constraint configurations are rejected", {
  dat <- simulate_re_mixture(60, c(0.5, 0.5), c(-1, 1), c(0.1, 0.1))
  expect_error(
    metamixture:::remix_cpp(dat$yi, dat$vi, K = 2L,
                           mu = "equal", tau2 = "equal"),
    "identical|unidentified"
  )
  expect_error(
    metamixture:::remix_cpp(dat$yi, dat$vi, K = 2L, mu = "bla"),
    "mu"
  )
  expect_error(
    metamixture:::remix_cpp(dat$yi, dat$vi, K = 2L, tau2 = "bla"),
    "tau2"
  )
})
