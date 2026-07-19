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
| `lassosir_fit()` | The Lasso-SIR comparison baseline — a thin wrapper around the external `LassoSIR` package (Lin, Zhao & Liu 2019), not a reimplementation. | Algorithm 1, Section 2.1 |
| `costeta()` | Cosine between two directions — the proximity criterion (PC). | Section 3.1 |
| `TPR.funct()`, `FPR.funct()`, `FDR.funct()` | Variable-selection accuracy metrics reported in Tables 1-6 (alongside `Cor`, computed directly with base R's `cor()`). | Section 3.1 |
| `X_rand()` | Predictor simulator: block-exchangeable (Setting 5), autoregressive (Settings 3-4), and a new spiked-eigenvalue option (Settings 1-2). | Section 3.1 |
| `normalize_vector()`, `normalize_matrix()` | Utility helpers. | — |

 

## About the two comparison baselines

The paper compares SIR-WLS against two baselines: plain Lasso, and
Lasso-SIR (Lin, Zhao & Liu 2019, *JASA*). Neither is reimplemented here —
Lasso-SIR in particular is someone else's published method with its own
CRAN package (`install.packages("LassoSIR")`). `lasso_fit()` and
`lassosir_fit()` are thin, consistently-named wrappers (matching exactly
the calls your scripts already use, e.g. `LassoSIR(X, y, H, no.dim = 1,
solution.path = FALSE, categorical = FALSE, nfolds = 10, screening =
FALSE)`) so all three methods compared in Tables 1-6 are callable from
`sirwls` with the same interface.

 

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
