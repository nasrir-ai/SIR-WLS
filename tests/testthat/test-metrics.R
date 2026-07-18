test_that("costeta returns 1 for parallel vectors and 0 for orthogonal ones", {
  expect_equal(costeta(c(1, 0), c(2, 0)), 1)
  expect_equal(costeta(c(1, 0), c(0, 5)), 0)
  expect_equal(costeta(c(1, 0), c(-3, 0)), -1)
})

test_that("selection metrics agree on a hand-worked example", {
  beta    <- c(1, 1, 0, 0, 0)
  hatbeta <- c(1, 0, 0.3, 0, 0)
  # true active set {1,2}; estimated active set {1,3}
  expect_equal(TPR.funct(beta, hatbeta), 0.5)     # TP = {1} -> 1/2
  expect_equal(FPR.funct(beta, hatbeta), 1 / 3)   # FP = {3} out of 3 true negatives
  expect_equal(FDR.funct(beta, hatbeta), 0.5)     # FP=1, TP=1 -> 1/(1+1)
})

test_that("selection metrics are perfect for an exact recovery", {
  beta <- c(1, 0, 2, 0)
  expect_equal(TPR.funct(beta, beta), 1)
  expect_equal(FPR.funct(beta, beta), 0)
  expect_equal(FDR.funct(beta, beta), 0)
})

test_that("selection metrics degrade gracefully for edge cases", {
  expect_equal(TPR.funct(rep(0, 4), rep(0, 4)), 0)  # no true positives to find
  expect_equal(FPR.funct(rep(1, 4), rep(1, 4)), 0)  # no true negatives to mislabel
  expect_equal(FDR.funct(c(1, 0, 0), c(0, 0, 0)), 0) # nothing selected -> FDR 0 by convention
})
