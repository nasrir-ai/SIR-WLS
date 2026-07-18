#' Internal: fetch dr:::dr.slices.arc regardless of export status
#'
#' The \pkg{dr} package's `dr.slices.arc()` is used unexported in some
#' package versions and exported in others. This helper fetches it from the
#' package namespace either way, so `wls.sir()` doesn't break across
#' \pkg{dr} versions.
#'
#' @param y Numeric response vector.
#' @param h Integer number of slices.
#' @return A list with `slice.indicator` and `slice.sizes`, as returned by
#'   `dr::dr.slices.arc()`.
#' @keywords internal
#' @noRd
.dr_slices_arc <- function(y, h) {
  if (!requireNamespace("dr", quietly = TRUE)) {
    stop("Package 'dr' is required for wls.sir(). Install it with install.packages(\"dr\").")
  }
  fn <- tryCatch(
    getExportedValue("dr", "dr.slices.arc"),
    error = function(e) get("dr.slices.arc", envir = asNamespace("dr"))
  )
  fn(y, h)
}

#' SIR-WLS: sliced inverse regression via weighted leverage score
#'
#' Algorithm 2 of Asrir and Mkhadri (2025). A model-free weighted leverage
#' score (Section 2.2, after Zhong et al. 2012) screens the columns of `x`
#' using the singular value decomposition of the (centered) design matrix;
#' classical sliced inverse regression is then applied only to the screened
#' predictors, and the resulting direction is extended with zeros for the
#' unselected predictors.
#'
#' @param x A numeric predictor matrix, `n` rows by `p` columns.
#' @param y A numeric response vector of length `n` (continuous), or a
#'   vector of class labels if `categorical = TRUE`.
#' @param nslice Integer number of slices used when `categorical = FALSE`
#'   (default 10); ignored when `categorical = TRUE`, where each distinct
#'   value of `y` is its own slice.
#' @param cn1 Numeric tuning constant for the (optional) BIC selection of
#'   the number of singular directions used to build the leverage scores,
#'   only used when `choose.dir = TRUE`.
#' @param cn2 Numeric tuning constant for the BIC selection of the number
#'   of screened predictors `p0` (Section 2.2, `A = {j : omega_j >=
#'   omega_(p0)}`).
#' @param choose.dir Logical. If `FALSE` (default, matching the published
#'   algorithm), all `r = min(n, p)` singular directions of `x` are used to
#'   build the weighted leverage scores. If `TRUE`, a BIC-type criterion
#'   first selects a smaller number of "spiked" directions to use instead
#'   — an extra option found in the author's exploratory code, not
#'   described in the published paper itself.
#' @param categorical Logical. Set to `TRUE` for a categorical/discrete `y`
#'   (each distinct value becomes its own slice), `FALSE` (default) for a
#'   continuous `y` sliced into `nslice` groups via `dr::dr.slices.arc()`.
#' @param ndim Integer number of EDR directions to estimate (default 1, the
#'   single-index case).
#'
#' @return A list with three components:
#' \describe{
#'   \item{wls}{Numeric vector of length `p`: the weighted leverage score
#'     `omega_j` for every predictor.}
#'   \item{select}{Integer vector: the indices of the `p0` predictors
#'     selected by the BIC criterion (the screened set `A`).}
#'   \item{betahat}{A `p x ndim` matrix: the estimated direction(s), zero
#'     outside `select`.}
#' }
#'
#' @details Predictors are screened by a two-step BIC procedure: first
#'   (only if `choose.dir = TRUE`) an eigenvalue-based BIC picks how many
#'   leading singular directions to use; then the per-predictor weighted
#'   leverage scores `omega_j` are ranked and a second BIC picks how many of
#'   the top-ranked predictors to keep (`p0`). Sliced inverse regression
#'   ([Rdimtools::do.rsir()] with `regmethod = "Ridge"`, for numerical
#'   stability on the screened subset) is then applied to `(x[, select], y)`
#'   to estimate the direction on the screened predictors, which is then
#'   zero-extended back to the original `p` dimensions.
#'
#' @references Asrir, N. and Mkhadri, A. (2025). A weighted leverage score
#'   approach versus Lasso sliced inverse regression: a comparative study.
#'   *Journal of Statistical Computation and Simulation*.
#'   \doi{10.1080/00949655.2025.2550340}. The weighted leverage score
#'   screening method is due to Zhong, W., Zhang, T., Zhu, Y., and Liu, J.
#'   S. (2012). Correlation pursuit: forward stepwise variable selection for
#'   index models. *Journal of the Royal Statistical Society Series B*,
#'   74(5), 849-870.
#'
#' @export
#' @examples
#' set.seed(1)
#' n <- 200; p <- 50
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- c(rep(1, 5), rep(0, p - 5))
#' y <- x %*% beta + rnorm(n)
#' fit <- wls.sir(x, y)
#' fit$select
#' costeta(beta, fit$betahat[, 1])
wls.sir <- function(x, y, nslice = 10, cn1 = 0.1, cn2 = 1, choose.dir = FALSE,
                     categorical = FALSE, ndim = 1) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  h <- nslice
  x <- scale(x, scale = FALSE)

  ## Slice y ----
  if (categorical) {
    index <- as.numeric(factor(y))
    nh <- as.numeric(summary(factor(y)))
    h <- length(nh)
  } else {
    slice <- .dr_slices_arc(y, h)
    index <- slice$slice.indicator
    nh <- slice$slice.sizes
  }

  ## SVD of the (centered) design matrix, r = min(n, p) directions ----
  svdx <- svd(x)
  u <- svdx$u
  d <- svdx$d
  v <- svdx$v

  if (!choose.dir) {
    dir <- min(n, p)
  } else {
    theta <- d^2 / (d[1])^2 + 1
    loglik <- penalty <- rep(0, length(d))
    for (i in seq_along(d)) {
      loglik[i] <- if (i < length(d)) {
        sum(log(theta[(i + 1):length(d)]) + 1 - theta[(i + 1):length(d)])
      } else 0
      penalty[i] <- i * cn1 / sqrt(n)
    }
    BIC_dir <- -loglik + penalty
    dir <- which.min(BIC_dir)
  }

  ## Weighted leverage score omega_j for every predictor ----
  w <- matrix(nrow = length(nh), ncol = dir)
  for (j in seq_len(dir)) {
    for (i in seq_along(nh)) w[i, j] <- sum(u[, j] * (index == i)) / nh[i]
  }

  uut <- array(dim = c(length(nh), dir, dir))
  for (i in seq_along(nh)) uut[i, , ] <- nh[i] * (w[i, ] %*% t(w[i, ]))
  sigma_slice <- colSums(uut)

  wls <- numeric(p)
  for (j in seq_len(p)) {
    vj <- v[j, seq_len(dir)]
    wls[j] <- as.numeric(t(vj) %*% sigma_slice %*% vj)
  }

  ## BIC selection of the number of screened predictors p0 ----
  wls_sort <- sort(wls, decreasing = TRUE)
  loglik <- penalty <- rep(0, min(n, p))
  for (k in seq_len(min(n, p))) {
    loglik[k] <- -log(sum(wls_sort[seq_len(k)]))
    penalty[k] <- (log(n) + cn2 * log(p)) * k / max(n, p)
  }
  BIC <- loglik + penalty
  sel_k <- which.min(BIC)
  select <- order(wls, decreasing = TRUE)[seq_len(sel_k)]

  ## SIR on the screened predictors, zero-extended back to p dimensions ----
  z <- x[, select, drop = FALSE]
  outsir <- Rdimtools::do.rsir(z, y, ndim = ndim, regmethod = "Ridge")

  betahat <- matrix(0, nrow = p, ncol = ndim)
  for (j in seq_len(ndim)) {
    betahat[select, j] <- Re(outsir$projection[, j])
  }

  list(wls = wls, select = select, betahat = betahat)
}
