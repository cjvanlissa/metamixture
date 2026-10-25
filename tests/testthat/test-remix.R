simulate_re_mixture <- function(n, pi, mu, tau2,
                                vi_range = c(0.005, 0.025)) {
  K <- length(pi)
  stopifnot(length(mu) == K, length(tau2) == K,
            all(pi > 0), abs(sum(pi) - 1) < 1e-10,
            all(tau2 >= 0))
  vi <- stats::runif(n, vi_range[1], vi_range[2])
  z <- sample.int(K, n, replace = TRUE, prob = pi)
  yi <- stats::rnorm(n, mean = mu[z], sd = sqrt(vi + tau2[z]))
  list(yi = yi, vi = vi, z = z)
}


testthat::test_that("one-class random-effects model recovers parameters", {
  set.seed(828)
  dat <- simulate_re_mixture(400, pi = 1, mu = 0.35, tau2 = 0.12)
  fit <- remix(dat$yi, dat$vi, classes = 1L)

  testthat::expect_length(OpenMx::omxGetParameters(fit, free = TRUE), 2L)
  testthat::expect_equal(fit$weights$values[1,1], 1, tolerance = 1e-12)
  testthat::expect_lt(abs(fit$M$values[1,] - 0.35), 0.01)
  testthat::expect_lt(abs(fit$S$values[1,] - 0.12), 0.01)
  testthat::expect_true(is.finite(as.numeric(fit$output$fit)))
})

testthat::test_that("two classes recover varying means and heterogeneity", {
  set.seed(831)
  dat <- simulate_re_mixture(600, pi = c(0.60, 0.40),
                             mu = c(-1.2, 1.4), tau2 = c(0.035, 0.16))
  fit <- remix(dat$yi, dat$vi, classes = 2L,
               mu = "varying", tau2 = "varying")

  testthat::expect_length(OpenMx::omxGetParameters(fit, free = TRUE), 5L)
  testthat::expect_equal(fit$M$values[1,], c(-1.2, 1.4), tolerance = 0.10)
  testthat::expect_true(all(abs(diag(fit$S$values) - c(0.035, 0.16)) < 0.010))
  pi <- metamixture:::remix_weights(fit)
  testthat::expect_true(all(abs(pi - c(0.60, 0.40)) < 0.01))
  testthat::expect_equal(sum(pi), 1, tolerance = 1e-10)
})

testthat::test_that("equal means and varying heterogeneity are recovered", {
  set.seed(14885)
  dat <- simulate_re_mixture(800, pi = c(0.60, 0.40),
                             mu = c(0.25, 0.25),
                             tau2 = c(0.015, 1.0),
                             vi_range = c(0.002, 0.015))
  fit <- remix(dat$yi, dat$vi, classes = 2L, mu = "equal",
               tau2 = "varying")


  testthat::expect_length(OpenMx::omxGetParameters(fit, free = TRUE), 4L)
  testthat::expect_equal(diff(fit$M$values[1,]), 0, tolerance = 1e-10)
  testthat::expect_lt(abs(fit$M$values[1,1] - 0.25), 0.02)
  testthat::expect_true(all(abs(sort(diag(fit$S$values)) - c(0.015, 1.0)) < 0.02))
  testthat::expect_equal(sort(metamixture:::remix_weights(fit)), c(.4, 0.60), tolerance = 0.05)
})

testthat::test_that("varying means and equal heterogeneity are recovered", {
  set.seed(953)
  dat <- simulate_re_mixture(600, pi = c(0.55, 0.45),
                             mu = c(-1.0, 1.2), tau2 = c(0.09, 0.09))
  fit <- remix(dat$yi, dat$vi, classes = 2L, mu = "varying",
               tau2 = "equal")
    testthat::expect_length(OpenMx::omxGetParameters(fit, free = TRUE), 4L)
  testthat::expect_equal(diff(diag(fit$S$values)), 0, tolerance = 1e-10)
  testthat::expect_true(all(abs(sort(fit$M$values[1,]) - c(-1, 1.2)) < 0.15))
  testthat::expect_true(all(abs(diag(fit$S$values) - 0.09) < 0.01))
  testthat::expect_true(all(abs(sort(metamixture:::remix_weights(fit)) - c(0.45, 0.55)) < 0.02))
})

testthat::test_that("unfitted OpenMx model has reproducible starting values", {
  set.seed(67)
  dat <- simulate_re_mixture(35, c(0.6, 0.4), c(-1, 1),
                             c(0.05, 0.10))
  set.seed(777)
  a <- remix(dat$yi, dat$vi, classes = 2L, run = FALSE)
  set.seed(777)
  b <- remix(dat$yi, dat$vi, classes = 2L, run = FALSE)
  testthat::expect_equal(OpenMx::omxGetParameters(a),
                         OpenMx::omxGetParameters(b))
})

testthat::test_that("invalid constraint configurations are rejected", {
  set.seed(123)
  dat <- simulate_re_mixture(60, c(0.5, 0.5), c(-1, 1), c(0.1, 0.1))
  testthat::expect_error(
    remix(dat$yi, dat$vi, classes = 2L,
          mu = "equal", tau2 = "equal", run = FALSE),
    regexp = "identical|unidentified|equal|identif"
  )
  testthat::expect_error(
    remix(dat$yi, dat$vi, classes = 2L, mu = "bla", run = FALSE),
    "arg|mu"
  )
  testthat::expect_error(
    remix(dat$yi, dat$vi, classes = 2L, tau2 = "bla", run = FALSE),
    "arg|tau2"
  )
})

