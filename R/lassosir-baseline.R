#' Lasso-SIR baseline (via the external \pkg{LassoSIR} package)
#'
#' A thin, consistently-named wrapper around `LassoSIR::LassoSIR()` (Lin,
#' Zhao, and Liu 2019), the second comparison baseline used throughout the
#' paper (Section 2.1, Algorithm 1). Lasso-SIR is not reimplemented here —
#' it is someone else's published method with its own CRAN package; this
#' wrapper just calls it with the arguments the paper's own scripts use, so
#' that all three compared methods (`lasso_fit()`, `lassosir_fit()`,
#' [wls.sir()]) share a common calling convention within `sirwls`.
#'
#' @param x A numeric predictor matrix, `n` rows by `p` columns.
#' @param y A numeric response vector of length `n` (continuous), or a
#'   vector of class labels if `categorical = TRUE`.
#' @param H Integer number of slices (default 10).
#' @param ndim Integer number of EDR directions to estimate (default 1,
#'   passed as `LassoSIR::LassoSIR()`'s `no.dim`).
#' @param categorical Logical, passed through to `LassoSIR::LassoSIR()`
#'   (default `FALSE`).
#' @param screening Logical, passed through to `LassoSIR::LassoSIR()`
#'   (default `FALSE`, matching the paper's own simulation scripts).
#' @param nfolds Integer number of cross-validation folds (default 10).
#'
#' @return A `p x ndim` numeric matrix: the estimated direction(s).
#'
#' @details Requires the \pkg{LassoSIR} package (CRAN); install with
#'   `install.packages("LassoSIR")`. Reproduces exactly the call used in the
#'   author's own comparison scripts: `LassoSIR(x, y, H, no.dim = ndim,
#'   solution.path = FALSE, categorical = categorical, nfolds = nfolds,
#'   screening = screening)`.
#'
#' @references Lin, Q., Zhao, Z., and Liu, J. S. (2019). Sparse sliced
#'   inverse regression via Lasso. *Journal of the American Statistical
#'   Association*, 114(528), 1726-1739. \doi{10.1080/01621459.2018.1520115}.
#'   Used as a comparison baseline in Asrir, N. and Mkhadri, A. (2025). A
#'   weighted leverage score approach versus Lasso sliced inverse
#'   regression: a comparative study. *Journal of Statistical Computation
#'   and Simulation*. \doi{10.1080/00949655.2025.2550340}.
#'
#' @seealso [wls.sir()], [lasso_fit()]
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("LassoSIR", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 300; p <- 50
#'   x <- matrix(rnorm(n * p), n, p)
#'   beta <- c(rep(1, 5), rep(0, p - 5))
#'   y <- x %*% beta + rnorm(n)
#'   beta.hat <- lassosir_fit(x, y)
#'   costeta(beta, beta.hat[, 1])
#' }
#' }
lassosir_fit <- function(x, y, H = 10, ndim = 1, categorical = FALSE,
                          screening = FALSE, nfolds = 10) {
  if (!requireNamespace("LassoSIR", quietly = TRUE)) {
    stop("Package 'LassoSIR' is required for lassosir_fit(). Install it with install.packages(\"LassoSIR\").")
  }
  fit <- LassoSIR::LassoSIR(
    x, y, H = H, no.dim = ndim, solution.path = FALSE,
    categorical = categorical, nfolds = nfolds, screening = screening
  )
  as.matrix(fit$beta)
}
