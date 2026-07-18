# sirwls

Companion R package for:

> Asrir, N. and Mkhadri, A. (2025). A weighted leverage score approach
> versus Lasso sliced inverse regression: a comparative study. *Journal of
> Statistical Computation and Simulation*.
> <https://doi.org/10.1080/00949655.2025.2550340>

## What this package contains

| Function | Purpose | Paper reference |
|---|---|---|
| `wls.sir()` | **SIR-WLS**, the main estimator: weighted-leverage-score screening (SVD-based), then SIR on the screened predictors, zero-filled elsewhere. | Algorithm 2, Section 2.3 |
| `lasso_fit()` | The plain Lasso comparison baseline. | Section 3, "usual Lasso" |
| `costeta()` | Cosine between two directions — the proximity criterion (PC). | Section 3.1 |
| `TPR.funct()`, `FPR.funct()`, `FDR.funct()` | Variable-selection accuracy metrics reported in Tables 1-6 (alongside `Cor`, computed directly with base R's `cor()`). | Section 3.1 |
| `X_rand()` | Predictor simulator: block-exchangeable (Setting 5), autoregressive (Settings 3-4), and a new spiked-eigenvalue option (Settings 1-2). | Section 3.1 |
| `normalize_vector()`, `normalize_matrix()` | Utility helpers. | — |

Unlike `nctsir` (the companion package for your other paper, where the
paper explicitly disclaims variable selection), this paper's whole point is
variable selection, and it reports TPR/FPR/FDR throughout Tables 1-6 — so
those three are included here as first-class, paper-verified metrics
(`FNR`/`TNR` are not reported anywhere in this paper, so they're left out,
same "only ship what's actually used" rule as before).

## Verified against the paper

- **`wls.sir()` matches Algorithm 2 and Section 2.2's weighted leverage
  score exactly**: SVD of the centered design matrix, per-predictor
  quadratic-form leverage scores from the slice-mean structure of the left
  singular vectors, a BIC criterion to pick how many predictors to keep
  (`p0`), SIR on the screened subset, zero-fill elsewhere.
- **Good news this time: no numerical bugs found in the core algorithm.**
  Both copies of `wls.sir()` in `myfunctions.R` were correct and identical.
  The only real bugs were in code *around* it — see below.
- **`X_rand(type = "spiked")` is new**: Settings 1-2 use a spiked-eigenvalue
  covariance (`Lambda = diag(81 or 51 decreasing spikes, 1, ..., 1)`,
  `x_i = Lambda^{1/2} u_i`) that I reconstructed directly from the paper's
  explicit formula in Section 3.1, since the author's own scripts only had
  an incomplete, commented-out attempt at it (missing the square root, and
  never wired into the main simulation loop).
- **The metrics match Section 3.1's definitions exactly**: PC (mean
  absolute cosine), Cor (correlation of the two projections — plain
  `cor()`, no wrapper needed), TPR, FPR, and FDR = FP / (FP + TP) —
  implemented directly from `beta`/`hatbeta`, rather than the
  `(p-q)*FPR/((p-q)*FPR+q*TPR)` approximation used ad hoc in some of the
  exploratory scripts (which needs `p`/`q` supplied and assumes a
  homogeneous split; the direct formula doesn't).

## Where the 5 uploaded files fit in

- **`myfunctions.R`** — the same giant personal file behind `nctsir`
  (~5,100 lines, several unrelated projects). `wls.sir()` (defined
  identically twice, lines ~3324 and ~3891) is the one this paper actually
  uses. Also found there but **excluded**:
    - `wls.sir.fdr()` — genuinely broken: references `t`, `Ta`, `alf`, `W`,
      `q`, none of which are defined inside the function; would crash if
      called. Not used by any of your three WLS scripts.
    - `cen.wls()` / `cen.wls.sir()` / `cen.wls.cox()` — a censored/survival
      data extension of WLS for a different project. `cen.wls()` is also
      duplicated and corrupted mid-file (a second definition starts before
      the first one's closing brace), but since none of it is used here it
      doesn't affect this package.
    - `AdaptCHOMPwithPIC()` ("SirCHOMP") — an unrelated method compared in
      early exploratory scripts, dropped from the final published
      comparison (Lasso / Lasso-SIR / SIR-WLS only).
- **`sir_wls_run_example.R`** — the script structurally closest to the
  published Tables 1-6: compares Lasso, Lasso-SIR, and SIR-WLS by PC, Cor,
  TPR, FPR, FDR, and even generates LaTeX table output. Kept as reference
  under `inst/scripts/sir_wls_run_example.R`.
- **`wls.sirvs nct.Sir vs lassoSIR .R`** — an earlier/parallel exploratory
  comparison that also includes `sir.nct` (your *other* paper's method) and
  `SirCHOMP`, neither of which appear in this paper's published comparison.
  Kept as reference (`inst/scripts/wls-sir-vs-nctsir-vs-lassoSIR.R`), not
  treated as the paper's own simulation code.
- **`real app wls .R`** — contains the paper's three real-data examples
  (Riboflavin, `~line 288`; CHIN, `~line 722`; SRBCT/Khan, `~line 1017`),
  plus two extra explorations not in the paper (Arcene, Wine datasets).
  Kept as reference (`inst/scripts/real-app-wls.R`); all four scripts use
  hard-coded local file paths (e.g. `~/Desktop/phd /R codes /...`) so none
  are directly runnable as-is.

## Still open

- I don't have R installed in this sandbox, so — same as with `nctsir` —
  I could not run `devtools::check()` or the test suite myself. Everything
  was checked by hand for syntax and cross-checked against the PDF.
- None of the three real datasets (Riboflavin via the `hdi` package, CHIN
  via `datamicroarray`, SRBCT via `plsgenomics`) are bundled — they're all
  available from CRAN/GitHub packages named in `real-app-wls.R`, but I
  haven't tried to build a clean, runnable vignette from that script since
  it's full of local paths and exploratory dead ends. Say the word if you
  want me to build one.
- The `choose.dir = TRUE` option in `wls.sir()` (an extra BIC step to pick
  a reduced number of singular directions before computing the leverage
  scores) is present in the source but not described in the published
  algorithm itself (which always uses `r = min(n, p)`, i.e.
  `choose.dir = FALSE`). It's kept as an option, documented as such, but
  not required for reproducing the paper.

## Verifying it yourself

```r
install.packages(c("devtools", "roxygen2", "dr", "Rdimtools", "glmnet", "mvtnorm"))
devtools::document("path/to/sirwls")
devtools::load_all("path/to/sirwls")
devtools::test("path/to/sirwls")
devtools::check("path/to/sirwls")
```

To cite the method, run `citation("sirwls")` once installed (see
`inst/CITATION`).
