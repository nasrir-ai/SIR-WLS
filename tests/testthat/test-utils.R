test_that("normalize_vector produces unit length and handles the zero vector", {
  v <- c(3, 4)
  nv <- normalize_vector(v)
  expect_equal(sqrt(sum(nv^2)), 1)
  expect_equal(normalize_vector(c(0, 0, 0)), c(0, 0, 0))
})

test_that("normalize_matrix normalizes columns and leaves an all-zero column alone", {
  mat <- matrix(c(3, 4, 0, 0, 1, 0), nrow = 2)
  out <- normalize_matrix(mat)
  expect_equal(sqrt(colSums(out^2)), c(1, 0, 1))
  expect_false(anyNA(out))
})
