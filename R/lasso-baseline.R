#' Plain Lasso baseline
#'
#' The "usual Lasso" comparison baseline used throughout the paper (Section
#' 3): 10-fold cross-validated lasso via \pkg{glmnet}, coefficients at
#' `lambda.min`, intercept dropped.
#'
#' @param x A numeric predictor matrix, `n` rows by `p` columns.
#' @param y A numeric response vector of length `n`.
#' @param family Passed to [glmnet::glmnet()]; `"gaussian"` (default) for a
#'   continuous response, `"binomial"` for a 0/1 response.
#' @param nfolds Integer number of cross-validation folds (default 10).
#'
#' @return A numeric vector of length `p`: the estimated coefficients
#'   (intercept excluded).
#'
#' @details This is a corrected, narrowed extraction of the author's
#'   original `Predict()` helper, which dispatched on a `method` string
#'   argument (`"Lasso"`, `"ALASSO-lasso"`, `"ALASSO-ridge"`,
#'   `"Elastic_net"`, `"SCAD"`, `"MCP"`) covering several unrelated
#'   projects. Two issues made the original risky to reuse as-is: the
#'   dispatch was on an exact string match against `"Lasso"`
#'   (title case) with no `match.arg()`/normalization, so a caller passing
#'   `method = "LASSO"` (as at least one of the author's own scripts did)
#'   silently fell through every branch and then errored on the undefined
#'   `fit` object; and the function had no default/error branch, so any
#'   unrecognised `method` failed the same way. Since this paper's Table 1-6
#'   "Lasso" column only ever uses the plain-lasso branch, this function
#'   keeps just that path under a fixed, unambiguous name instead of
#'   reproducing the dispatch bug.
#'
#' @seealso [wls.sir()]
#' @export
#' @examples
#' set.seed(1)
#' n <- 200; p <- 50
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- c(rep(1, 5), rep(0, p - 5))
#' y <- x %*% beta + rnorm(n)
#' beta.hat <- lasso_fit(x, y)
#' costeta(beta, beta.hat)
lasso_fit <- function(x, y, family = "gaussian", nfolds = 10) {
  x <- as.matrix(x)
  cvfit <- glmnet::cv.glmnet(x, y, alpha = 1, family = family, nfolds = nfolds)
  fit <- glmnet::glmnet(x, y, alpha = 1, family = family, lambda = cvfit$lambda.min)
  beta_coef <- stats::coef(fit, exact = TRUE)
  as.numeric(beta_coef[-1])
}
