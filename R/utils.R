#' Normalize a vector to unit length
#'
#' @param v A numeric vector.
#'
#' @return `v` rescaled to unit Euclidean length; `v` itself, unchanged, if
#'   it is exactly the zero vector.
#'
#' @export
#' @examples
#' normalize_vector(c(3, 4))
normalize_vector <- function(v) {
  magnitude <- sqrt(sum(v^2))
  if (magnitude == 0) return(v)
  v / magnitude
}

#' Normalize each column of a matrix to unit length
#'
#' @param mat A numeric matrix.
#'
#' @return A matrix of the same dimensions as `mat`, with each column
#'   rescaled to unit Euclidean length by [normalize_vector()].
#'
#' @export
#' @examples
#' normalize_matrix(matrix(1:6, nrow = 2))
normalize_matrix <- function(mat) {
  apply(mat, 2, normalize_vector)
}
