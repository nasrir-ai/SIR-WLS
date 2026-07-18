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
