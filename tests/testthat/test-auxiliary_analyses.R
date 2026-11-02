test_that("auxiliary methods work with remix models", {

  skip_if_not_installed("OpenMx")
  skip_if_not_installed("metamixture")
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("mice")

  library(tidySEM)
  set.seed(1234)

  # Simulate two classes
  n <- 150
  C <- sample(1:2, n, replace = TRUE)

  vi <- runif(n, 0.01, 0.05)
  yi <- rnorm(
    n,
    mean = c(-1.5, 1.5)[C],
    sd = sqrt(vi + 0.1)
  )

  # Distal outcome
  dat <- data.frame(
    y = rnorm(n, mean = c(0, 2)[C], sd = 1)
  )

  # Fit remix
  fit <- metamixture::remix(
    yi = yi,
    vi = vi,
    classes = 2,
    mu = "varying",
    tau2 = "equal"
  )

  # Pseudo-class
  pc <- tidySEM::pseudo_class(
    x = fit,
    model = lm(y ~ class, data = data),
    data = dat,
    m = 5
  )
  expect_equal(as.character(pc$term), c("(Intercept)", "class"))


  # BCH
  bch <- tidySEM::BCH(
    x = fit,
    model = "y ~ 1",
    data = dat
  )

  expect_s4_class(bch, "MxModel")
  expect_true(is.finite(bch@output$fit))
  expect_equal(sum(grepl("class\\d\\.M", names(coef(bch)))), 2)

  # Two-step maximum likelihood
  ts <- tidySEM::two_step(
    x = fit,
    model = "y ~ 1",
    data = dat
  )

  expect_s4_class(ts, "MxModel")
  expect_true(is.finite(ts@output$fit))
  expect_equal(sum(grepl("aux\\d_M", names(coef(ts)))), 2)

  # Corrected standard errors
  se <- attr(ts, "se_robust")

  expect_true(is.numeric(se))
  expect_true(length(se) > 0L)
  expect_true(all(is.finite(se)))
  expect_true(all(se >= 0))
})
