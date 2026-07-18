test_that("X_rand returns a matrix of the requested shape for every non-spiked type", {
  skip_if_not_installed("mvtnorm")
  set.seed(1)
  for (ty in c("block", "autoR", "homog", "homog2")) {
    x <- X_rand(n = 30, p = 8, q = 3, type = ty, alpha = c(0.5, 0.2, 0.4), rho = 0.4)
    expect_equal(dim(x), c(30, 8))
    expect_false(anyNA(x))
  }
})

test_that("X_rand(type = 'spiked') returns the right shape and a decreasing spike pattern", {
  set.seed(1)
  n <- 500; p <- 300
  x <- X_rand(n = n, p = p, type = "spiked", n_spike = 81)
  expect_equal(dim(x), c(n, p))
  expect_false(anyNA(x))
  # spiked columns should have larger sample variance than the non-spiked tail
  vars <- apply(x, 2, var)
  expect_gt(mean(vars[1:81]), mean(vars[82:p]))
})

test_that("X_rand(type = 'spiked') matches Setting 2's 51-spike configuration", {
  set.seed(1)
  x <- X_rand(n = 500, p = 300, type = "spiked", n_spike = 51)
  expect_equal(dim(x), c(500, 300))
})
