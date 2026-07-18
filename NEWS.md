# sirwls 0.1.0

* Initial package extracted from the author's working files
  (`myfunctions.R`, shared with the `nctsir` package, plus three
  WLS-specific analysis scripts) for the paper "A weighted leverage score
  approach versus Lasso sliced inverse regression: a comparative study".
* Exported: `wls.sir()`, `lasso_fit()`, `costeta()`, `TPR.funct()`,
  `FPR.funct()`, `FDR.funct()`, `normalize_vector()`, `normalize_matrix()`,
  `X_rand()` (with a new `type = "spiked"` covariance option reconstructed
  directly from the paper's Settings 1-2 formula, which was present only as
  incomplete/commented-out exploratory code in the author's scripts).
* `wls.sir()`'s core computation (weighted leverage score screening + SIR
  on the screened subset) was already correct in the source file — no
  numerical bugs found there, unlike the author's other two packages.
* Added a defensive lookup for `dr:::dr.slices.arc()` (`.dr_slices_arc()`),
  since that function is unexported in some versions of the `dr` package.
* Extracted `lasso_fit()` from the author's `Predict()` dispatcher, which
  covered several unrelated methods (`ALASSO-lasso`, `ALASSO-ridge`,
  `Elastic_net`, `SCAD`, `MCP`) behind a `method` string with no
  `match.arg()`/case-normalization and no error branch — at least one of
  the author's own scripts called it with `method = "LASSO"` (wrong case),
  which silently fell through every branch and crashed on the undefined
  `fit` object. Since only the plain-Lasso branch is used by this paper,
  `lasso_fit()` keeps just that path, unambiguously, correctly cased.
* Excluded `wls.sir.fdr()` (unused elsewhere, and genuinely broken: it
  references `t`, `Ta`, `alf`, `W`, `q` that are never defined in that
  function and would error if called) and the `cen.wls*()` family (a
  censored/survival-data extension of WLS for a different, unrelated
  project — also found duplicated and corrupted mid-file in
  `myfunctions.R`, an unrelated issue that doesn't affect this package
  since none of that code was used here).
* Excluded `AdaptCHOMPwithPIC()`/"SirCHOMP" (an unrelated method the
  author explored for comparison in early drafts, not part of the
  published paper's Lasso vs. Lasso-SIR vs. SIR-WLS comparison).
