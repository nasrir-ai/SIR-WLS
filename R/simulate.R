#' Simulate correlated Gaussian predictors
#'
#' Generates an `n x p` matrix of predictors under one of the covariance
#' structures used in the paper's five simulation settings (Section 3.1).
#'
#' @param n Integer sample size.
#' @param p Integer number of predictors.
#' @param q Integer number of "active" predictors (used by `type = "block"`
#'   and `type = "homog2"` only — irrelevant to `type = "spiked"`, where the
#'   active set is chosen independently of the covariance structure; see
#'   Details).
#' @param type One of `"block"` (Setting 5), `"autoR"` (Settings 3-4,
#'   `Cov(x_i, x_j) = rho^|i-j|`), `"spiked"` (Settings 1-2), `"homog"`,
#'   `"homog2"`.
#' @param alpha Length-3 vector \eqn{(\alpha_1, \alpha_2, \alpha_3)} used
#'   when `type = "block"`: within-active-block, active/inactive
#'   cross-block, and within-inactive-block correlation, respectively —
#'   exactly Section 3.1's block-exchangeable \eqn{\Sigma_{11}},
#'   \eqn{\Sigma_{12}}, \eqn{\Sigma_{22}}.
#' @param rho Correlation parameter used by `type` `"autoR"`, `"homog"`,
#'   `"homog2"`.
#' @param S Integer indices (length `q`) of the active predictors, used by
#'   `type = "homog2"`. Defaults to the first `q` predictors.
#' @param n_spike Integer, used only by `type = "spiked"`: the number of
#'   spiked eigenvalues (81 for Setting 1, 51 for Setting 2 of the paper).
#'
#' @return An `n x p` numeric matrix.
#'
#' @details
#' For `type = "spiked"`, predictors are generated as \eqn{x_i =
#' \Lambda^{1/2} u_i}, \eqn{u_i \sim N(0, I_p)}, where \eqn{\Lambda =
#' \mathrm{diag}(\lambda_{spike}, \ldots, \lambda_{spike} - n\_spike + 1,
#' 1, \ldots, 1)} with `lambda_spike = n_spike - 1 + ceiling(p / sqrt(n))`,
#' matching Settings 1-2 of Section 3.1 exactly (`n_spike = 81` and `51`
#' respectively). Unlike the other `type`s, the active-predictor index set
#' for `"spiked"` (`S = {1, 10, 15, 20, 25, 30}` in Setting 1, extended with
#' `{40, 45, 50}` in Setting 2) is unrelated to the covariance structure —
#' build `Beta` directly with those indices, the same way the paper does.
#'
#' @references Section 3.1 in Asrir, N. and Mkhadri, A. (2025). A weighted
#'   leverage score approach versus Lasso sliced inverse regression: a
#'   comparative study. *Journal of Statistical Computation and
#'   Simulation*. \doi{10.1080/00949655.2025.2550340}.
#'
#' @export
#' @examples
#' x <- X_rand(n = 50, p = 10, q = 3, type = "block", alpha = c(0.5, 0.2, 0.4))
#' dim(x)
#'
#' # Setting 1 of the paper (81 spiked eigenvalues)
#' x1 <- X_rand(n = 500, p = 300, type = "spiked", n_spike = 81)
#' s <- c(1, 10, 15, 20, 25, 30)
#' Beta <- rep(0, 300); Beta[s] <- 1
X_rand <- function(n, p, q = NULL, type = c("block", "autoR", "spiked", "homog", "homog2"),
                    alpha = c(0.5, 0.2, 0.4), rho = 0, S = if (!is.null(q)) seq_len(q) else NULL,
                    n_spike = 81) {
  type <- match.arg(type)

  if (type == "spiked") {
    lambda_top <- (n_spike - 1) + ceiling(p / sqrt(n))
    spikes <- lambda_top - (seq_len(n_spike) - 1)
    eigvals <- c(spikes, rep(1, p - n_spike))
    u <- matrix(stats::rnorm(n * p), n, p)
    return(u %*% diag(sqrt(eigvals)))
  }

  sigma <- switch(
    type,
    block = {
      C11 <- matrix(alpha[1], q, q); diag(C11) <- 1
      C22 <- matrix(alpha[3], p - q, p - q); diag(C22) <- 1
      C12 <- matrix(alpha[2], q, p - q)
      rbind(cbind(C11, C12), cbind(t(C12), C22))
    },
    autoR = rho^abs(outer(seq_len(p), seq_len(p), "-")),
    homog = {
      m <- matrix(rho, p, p)
      diag(m) <- 1
      m
    },
    homog2 = {
      m <- matrix(0, p, p)
      m[S, S] <- rho
      m[-S, -S] <- rho
      m[S, -S] <- 0.1
      m[-S, S] <- 0.1
      diag(m) <- 1
      m
    }
  )

  if (!isSymmetric(sigma)) sigma <- 0.5 * (sigma + t(sigma))

  mvtnorm::rmvnorm(n, rep(0, p), sigma)
}
