test_that("cv_enumeration selects classes using cross-validation", {

  set.seed(123)

  # Simulate meta-analytic data
  dat <- sim_remix_fixed_tau(k = 100, es = 2)

  # Estimate candidate models
  models <- remix(
      yi = dat$yi,
      vi = dat$vi,
      classes = 1:3,
      tau2 = "equal"
    )

  # Five-fold cross-validation
  set.seed(456)
  cv <- cv_enumeration(models, folds = 5)

  # Check return structure
  expect_type(cv, "list")
  expect_named(
    cv,
    c("best_class", "scores", "fold_scores", "fold_id")
  )

  # Selected model is one of the candidates
  expect_length(cv$best_class, 1L)
  expect_true(cv$best_class %in% seq_along(models))

  # Predictive scores
  expect_length(cv$scores, length(models))
  expect_true(all(is.finite(cv$scores)))

  # Fold-level scores
  expect_true(is.matrix(cv$fold_scores))
  expect_equal(dim(cv$fold_scores), c(5L, 3L))
  expect_true(all(is.finite(cv$fold_scores)))

  # Total scores equal sum of fold scores
  expect_equal(
    unname(cv$scores),
    unname(colSums(cv$fold_scores))
  )

  # Best model has lowest predictive score
  expect_equal(
    unname(cv$best_class),
    unname(which.min(cv$scores))
  )

  # Every observation belongs to exactly one fold
  expect_length(cv$fold_id, nrow(dat))
  expect_setequal(unique(cv$fold_id), seq_len(5))

  # Folds are approximately equally sized
  expect_lte(
    diff(range(tabulate(cv$fold_id))),
    1L
  )
})
