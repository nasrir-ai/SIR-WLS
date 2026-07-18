test_that("lasso_fit returns p coefficients and recovers a sparse signal reasonably well", {
  skip_if_not_installed("glmnet")
  set.seed(1)
  n <- 200; p <- 30
  x <- matrix(rnorm(n * p), n, p)
  beta <- c(rep(1, 4), rep(0, p - 4))
  y <- x %*% beta + rnorm(n, sd = 0.3)

  beta_hat <- lasso_fit(x, y)
  expect_length(beta_hat, p)
  expect_gt(abs(costeta(beta, beta_hat)), 0.8)
})

test_that("lasso_fit works for a binomial response", {
  skip_if_not_installed("glmnet")
  set.seed(1)
  n <- 200; p <- 20
  x <- matrix(rnorm(n * p), n, p)
  beta <- c(rep(1, 3), rep(0, p - 3))
  y <- rbinom(n, 1, plogis(x %*% beta))

  expect_no_error(beta_hat <- lasso_fit(x, y, family = "binomial"))
  expect_length(beta_hat, p)
})
