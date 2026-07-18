#' Cosine between two direction vectors (proximity criterion)
#'
#' The building block of the paper's "proximity criterion" (PC, Section
#' 3.1): \eqn{PC = \frac{1}{N}\sum_{k=1}^N |\cos(\hat\beta_k, \beta)|}, the
#' mean absolute cosine between the true direction and its estimate over
#' `N` replications.
#'
#' @param x,y Numeric vectors (or single-column matrices) of the same length.
#'
#' @return The cosine similarity between `x` and `y`, a value in `[-1, 1]`.
#'   Since sufficient-dimension-reduction directions are only identified up
#'   to sign, the paper reports `abs(costeta(beta, hatbeta))` averaged over
#'   replications.
#'
#' @references Section 3.1 in Asrir, N. and Mkhadri, A. (2025). A weighted
#'   leverage score approach versus Lasso sliced inverse regression: a
#'   comparative study. *Journal of Statistical Computation and
#'   Simulation*. \doi{10.1080/00949655.2025.2550340}.
#'
#' @export
#' @examples
#' costeta(c(1, 0), c(1, 1))
costeta <- function(x, y) {
  x <- as.matrix(x); y <- as.matrix(y)
  as.numeric(crossprod(x, y) / (norm(x, "f") * norm(y, "f")))
}

#' Variable-selection accuracy metrics (TPR, FPR, FDR)
#'
#' `TPR.funct`, `FPR.funct`, and `FDR.funct` compare the support (nonzero
#' pattern) of an estimate `hatbeta` against the true coefficient vector
#' `beta`. These are exactly the three variable-selection metrics reported
#' in the paper's Tables 1-6, alongside the proximity criterion
#' ([costeta()]) and the correlation `cor(X \%*\% beta, X \%*\% hatbeta)`
#' (for which base R's `cor()` is used directly — no wrapper needed).
#'
#' @param beta Numeric vector: the true coefficient vector.
#' @param hatbeta Numeric vector, same length as `beta`: the estimated
#'   coefficient vector.
#'
#' @return A single numeric value in `[0, 1]`.
#'
#' @details `TPR` is the proportion of truly active predictors that are
#'   selected; `FPR` is the proportion of truly inactive predictors that
#'   are (falsely) selected; `FDR` is the proportion of *selected*
#'   predictors that are false positives (`FP / (FP + TP)`), matching the
#'   paper's definition exactly ("the ratio of false positives to the total
#'   number of positive predictions", Section 3.1). This differs from the
#'   `(p - q) * FPR / ((p - q) * FPR + q * TPR)` formula used ad hoc in some
#'   of the author's exploratory scripts, which back-computes an
#'   *approximate* FDR from FPR/TPR assuming a homogeneous active/inactive
#'   split; `FDR.funct()` instead computes it directly from `beta` and
#'   `hatbeta`, which is exact and does not need `p`/`q` supplied
#'   separately. `FDR.funct()` returns 0 when nothing is selected (by
#'   convention, no false discoveries among zero discoveries).
#'
#' @references Section 3.1 in Asrir, N. and Mkhadri, A. (2025). A weighted
#'   leverage score approach versus Lasso sliced inverse regression: a
#'   comparative study. *Journal of Statistical Computation and
#'   Simulation*. \doi{10.1080/00949655.2025.2550340}.
#'
#' @name selection-metrics
#' @examples
#' beta <- c(1, 1, 0, 0, 0)
#' hatbeta <- c(1, 0, 0.3, 0, 0)
#' TPR.funct(beta, hatbeta)
#' FPR.funct(beta, hatbeta)
#' FDR.funct(beta, hatbeta)
NULL

#' @rdname selection-metrics
#' @export
TPR.funct <- function(beta, hatbeta) {
  q <- sum(beta != 0)
  if (q == 0) return(0)
  TP <- sum(beta != 0 & hatbeta != 0)
  TP / q
}

#' @rdname selection-metrics
#' @export
FPR.funct <- function(beta, hatbeta) {
  p <- length(beta)
  q <- sum(beta != 0)
  if ((p - q) == 0) return(0)
  FP <- sum(hatbeta != 0 & beta == 0)
  FP / (p - q)
}

#' @rdname selection-metrics
#' @export
FDR.funct <- function(beta, hatbeta) {
  TP <- sum(beta != 0 & hatbeta != 0)
  FP <- sum(hatbeta != 0 & beta == 0)
  if ((TP + FP) == 0) return(0)
  FP / (TP + FP)
}
