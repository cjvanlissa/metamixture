test_that("num_classes works for latent profile models", {

  skip_if_not_installed("OpenMx")

  set.seed(456)

  n <- 100
  vi <- runif(n, 0.01, 0.05)
  yi <- rnorm(n)


  for (K in 1:2) {

  # Unestimated model
    model <- suppressWarnings(metamixture::remix(
      yi = yi,
      vi = vi,
      classes = K,
      run = FALSE
    ))

    expect_equal(metamixture:::num_classes(model), K)

    # Estimated model
    fit <- suppressWarnings(OpenMx::mxRun(
      model,
      intervals = FALSE,
      silent = TRUE
    ))

    expect_equal(metamixture:::num_classes(fit), K)
  }
})


test_that("num_classes works for lists of mixture models", {


  skip_if_not_installed("OpenMx")

  set.seed(456)

  n <- 100
  vi <- runif(n, 0.01, 0.05)
  yi <- rnorm(n)


  model <- suppressWarnings(metamixture::remix(
    yi = yi,
    vi = vi,
    classes = 1:2,
    run = FALSE
  ))
  expect_equal(unname(sapply(model, metamixture:::num_classes, USE.NAMES = F)), c(1L:2L))

  # Estimated model
  fit <- suppressWarnings(lapply(model, OpenMx::mxRun, intervals = FALSE, silent = TRUE))

  expect_equal(unname(sapply(fit, metamixture:::num_classes, USE.NAMES = F)), c(1L:2L))


})
