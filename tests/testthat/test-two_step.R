test_that("two_step auxiliary likelihood uses probability units", {

  skip_if_not_installed("OpenMx")

  set.seed(123)

  n <- 100
  z <- rnorm(n, mean = 1, sd = 2)

  aux <- as_ram(
    "z ~ 1",
    data = data.frame(z = z)
  )

  aux <- OpenMx::mxModel(
    aux,
    OpenMx::mxFitFunctionML(vector = TRUE)
  )

  aux <- suppressWarnings(OpenMx::mxRun(
    aux,
    silent = TRUE
  ))

  likelihood <- as.numeric(
    OpenMx::mxEvalByName(
      "fitfunction",
      aux
    )
  )

  manual <- dnorm(
    z,
    mean = aux$M$values[1, "z"],
    sd = sqrt(aux$S$values["z", "z"])
  )

  expect_equal(
    likelihood,
    manual,
    tolerance = 1e-6
  )
})
