test_that("sandwich_remix is invariant to exact duplication within clusters", {
  set.seed(1)

  n <- 100
  df <- metamixture::sim_remix_fixed_tau(n, es = 2)

  # Original sample: independent observations
  res1 <- metamixture::remix(
    yi = df$yi,
    vi = df$vi,
    classes = 2,
    tau2 = "equal"
  )

  sw1 <- metamixture::sandwich_remix(
    res1,
    cluster = seq_len(n),
    type = "CR0"
  )

  # Duplicate observations exactly, retaining cluster membership
  df2 <- rbind(df, df)
  df2$study <- rep(seq_len(n), 2)

  res2 <- metamixture::remix(
    yi = df2$yi,
    vi = df2$vi,
    classes = 2,
    tau2 = "equal"
  )

  sw2 <- metamixture::sandwich_remix(
    res2,
    cluster = df2$study,
    type = "CR0"
  )

  # Extract cluster-robust covariance matrices
  V1 <- attr(sw1, "sandwich")$vcov
  V2 <- attr(sw2, "sandwich")$vcov

  expect_true(is.matrix(V1))
  expect_true(is.matrix(V2))
  expect_true(all(is.finite(V1)))
  expect_true(all(is.finite(V2)))

  # Exact duplication should preserve the CR0 covariance
  expect_equal(
    V1,
    V2,
    tolerance = 1e-3
  )
  res2se <- res2@output$standardErrors
  tabres <- table_results(sw2, format_numeric = FALSE, columns = NULL)
  tabres <- tabres[match(rownames(res2se), tabres$openmx_label), ]
  expect_true(all(tabres$se - res2se >= 0))
})
