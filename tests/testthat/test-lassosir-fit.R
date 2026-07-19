test_that("lassosir_fit returns a p x ndim matrix and recovers a sparse signal reasonably well", {
  skip_if_not_installed("LassoSIR")
  set.seed(1)
  n <- 300; p <- 50
  x <- matrix(rnorm(n * p), n, p)
  beta <- c(rep(1, 5), rep(0, p - 5))
  y <- x %*% beta + rnorm(n)

  beta_hat <- lassosir_fit(x, y)
  expect_equal(dim(beta_hat), c(p, 1))
  expect_gt(abs(costeta(beta, beta_hat[, 1])), 0.5)
})

test_that("lassosir_fit errors clearly when LassoSIR is not installed", {
  skip_if(requireNamespace("LassoSIR", quietly = TRUE), "LassoSIR is installed, can't test the missing-package path")
  expect_error(lassosir_fit(matrix(1, 2, 2), c(1, 2)), "LassoSIR")
})
