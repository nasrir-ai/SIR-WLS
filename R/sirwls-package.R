#' sirwls: Sliced Inverse Regression via Weighted Leverage Score
#'
#' @description
#' Companion package for Asrir and Mkhadri (2025), "A weighted leverage
#' score approach versus Lasso sliced inverse regression: a comparative
#' study", *Journal of Statistical Computation and Simulation*,
#' \doi{10.1080/00949655.2025.2550340}.
#'
#' Provides the SIR-WLS estimator ([wls.sir()], Algorithm 2 of the paper), a
#' plain Lasso baseline ([lasso_fit()]), the proximity criterion
#' ([costeta()]) and variable-selection accuracy metrics ([TPR.funct()],
#' [FPR.funct()], [FDR.funct()]) used throughout the paper's empirical
#' study, and a correlated Gaussian predictor simulator matching its five
#' simulation settings ([X_rand()], Section 3.1).
#'
#' @section Where these functions came from:
#' This package was assembled from the same personal working file
#' (`myfunctions.R`) used for the author's other SIR package, `nctsir`,
#' which accumulates code for several different, unrelated SIR projects.
#' Only the functions actually exercised by this paper's own analysis
#' scripts were extracted here — see `NEWS.md` and the package README for
#' full provenance notes and the bugs checked/fixed along the way.
#'
#' @keywords internal
"_PACKAGE"
