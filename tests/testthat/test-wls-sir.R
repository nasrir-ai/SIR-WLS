test_that("wls.sir recovers a sparse direction reasonably well and returns correct shapes", {
  skip_if_not_installed("dr")
  skip_if_not_installed("Rdimtools")
  set.seed(42)
  n <- 300; p <- 50
  x <- matrix(rnorm(n * p), n, p)
  beta <- c(rep(1, 5), rep(0, p - 5))
  beta <- beta / sqrt(sum(beta^2))
  y <- x %*% beta + rnorm(n, sd = 0.5)

  fit <- wls.sir(x, y)
  expect_length(fit$wls, p)
  expect_true(all(fit$select >= 1 & fit$select <= p))
  expect_equal(dim(fit$betahat), c(p, 1))
  expect_true(all(which(fit$betahat[, 1] != 0) %in% fit$select))
  # up to sign, the estimated direction should be reasonably aligned with beta
  expect_gt(abs(costeta(beta, fit$betahat[, 1])), 0.6)
})

test_that("wls.sir works with categorical y (binary response)", {
  skip_if_not_installed("dr")
  skip_if_not_installed("Rdimtools")
  set.seed(1)
  n <- 200; p <- 30
  x <- matrix(rnorm(n * p), n, p)
  beta <- c(rep(1, 4), rep(0, p - 4))
  y <- rbinom(n, 1, plogis(x %*% beta))

  expect_no_error(fit <- wls.sir(x, y, categorical = TRUE))
})

test_that("wls.sir's dr.slices.arc lookup does not error regardless of export status", {
  skip_if_not_installed("dr")
  set.seed(1)
  y <- rnorm(50)
  out <- sirwls:::.dr_slices_arc(y, 10)
  expect_true(!is.null(out$slice.indicator))
  expect_true(!is.null(out$slice.sizes))
})
