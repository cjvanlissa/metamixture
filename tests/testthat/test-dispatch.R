test_that("tidySEM dispatches remix methods", {
  set.seed(1)
  model <- remix(
    yi = c(-1, -0.5, 0.5, 1),
    vi = rep(0.1, 4),
    classes = 2,
    run = FALSE
  )

  expect_identical(attr(model, "tidySEM"), "remix")

  result <- tidySEM::random_starts(
    model,
    nstarts = 2,
    run = FALSE
  )

  expect_type(result, "list")
  expect_length(result, 2)
})
