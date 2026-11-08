test_that("one_step jointly estimates remix and auxiliary parameters", {

  set.seed(6883)

  n <- 100

  df <- sim_remix_fixed_tau(
    k = n,
    es = 2
  )

  # Auxiliary outcome related to true class
  df$z <- rnorm(
    n,
    mean = df$C - 1,
    sd = 1
  )

  # Fit measurement-only mixture
  fit <- remix(
    yi = df$yi,
    vi = df$vi,
    classes = 2,
    tau2 = "equal"
  )

  # One-step joint estimation
  res <- tidySEM::one_step(
    x = fit,
    data = df$z
  )

  expect_s4_class(res, "MxModel")
  expect_true(isTRUE(res@.wasRun))
  expect_true(is.finite(res@output$fit))

  # Measurement parameters remain free
  measurement_labels <- names(
    OpenMx::omxGetParameters(
      fit,
      free = TRUE
    )
  )

  joint_labels <- names(
    OpenMx::omxGetParameters(
      res,
      free = TRUE
    )
  )

  expect_true(all(
    measurement_labels %in% joint_labels
  ))

  # Auxiliary parameters are also free
  expect_true(any(grepl("^aux", joint_labels)))

  # No two-step correction is required
  expect_null(attr(res, "se_robust"))

  res_two <- tidySEM::two_step(
    x = fit,
    data = df$z
  )

  res_one <- tidySEM::one_step(
    x = fit,
    data = df$z
  )
  # Free measurement model gives better fit:
  expect_true(res_one@output$Minus2LogLikelihood < res_two@output$Minus2LogLikelihood)
})
