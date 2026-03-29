# nmfkc: Non-negative Matrix Factorization with Kernel Covariates

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![GitHub version](https://img.shields.io/github/r-package/v/ksatohds/nmfkc)](https://github.com/ksatohds/nmfkc)

**nmfkc** is an R package that extends Non-negative Matrix Factorization (NMF) by incorporating **covariates** using kernel methods. It supports advanced features like rank selection via cross-validation, time-series modeling (NMF-VAR), supervised classification (NMF-LAB), **structural equation modeling with equilibrium interpretation (NMF-SEM)**, and **mixed-effects modeling with random effects (NMF-RE)**.

# Installation

```r
# install.packages("remotes")
remotes::install_github("ksatohds/nmfkc")
library(nmfkc)
```

# Help and Usage

```r
browseVignettes("nmfkc")
ls("package:nmfkc")
?nmfkc
```

# Citation

```r
citation("nmfkc")
```

# Quick Example

```r
library(nmfkc)

# Decompose a matrix Y into basis X and coefficient B with rank = 2
X_true <- cbind(c(1, 0, 1), c(0, 1, 0))
B_true <- cbind(c(1, 0), c(0, 1), c(1, 1))
Y <- X_true %*% B_true

res <- nmfkc(Y, rank = 2, epsilon = 1e-6)
plot(res)     # Convergence plot
summary(res)  # Summary statistics
```

See `browseVignettes("nmfkc")` for detailed examples covering rank selection, kernel NMF, time-series, classification, NMF-SEM, and NMF-RE.

# Comparison with Standard NMF

| Feature | Standard NMF | nmfkc |
| :--- | :--- | :--- |
| **Handles covariates** | No | **Yes** (Linear / Kernel) |
| **Structural equation modeling** | No | **Yes** (NMF-SEM) |
| **Mixed-effects / Random effects** | No | **Yes** (NMF-RE) |
| **Classification** | No | **Yes** (NMF-LAB) |
| **Time series modeling** | No | **Yes** (NMF-VAR) |
| **Nonlinearity** | No | **Yes** (Kernel) |
| **Clustering support** | Limited | **Yes** (Hard/Soft) |
| **Rank selection / CV** | Limited (ad hoc) | **Yes** (Element-wise CV, Column-wise CV) |

# Statistical Model

The **nmfkc** package builds upon the standard NMF framework by incorporating external information (covariates):

$$Y(P,N) \approx X(P,Q) \times C(Q,R) \times A(R,N)$$

- **$Y$**: Observation matrix ($P$ features × $N$ samples)
- **$A$**: Covariate matrix ($R$ covariates × $N$ samples); defaults to identity (standard NMF)
- **$X$**: Basis matrix — learned latent patterns
- **$C$**: Parameter matrix — links covariates to latent structure

### Extensions

- **NMF-RE**: Adds unit-specific random effects $U$: $Y = X(\Theta A + U) + \mathcal{E}$, estimated via ridge-type BLUP with wild bootstrap inference.
- **NMF-SEM**: Models endogenous feedback: $Y_1 \approx X(\Theta_1 Y_1 + \Theta_2 Y_2)$, with equilibrium mapping $(I - X\Theta_1)^{-1} X\Theta_2$.

# Main Functions

| Function | Description |
|:---|:---|
| `nmfkc()` | Core NMF with covariates ($Y \approx XCA$); supports kernel matrices and formula interface |
| `nmfre()` / `nmfre.inference()` | NMF with Random Effects + wild bootstrap inference |
| `nmf.sem()` / `nmf.sem.inference()` | NMF Structural Equation Model + inference for path coefficients |
| `nmfae()` / `nmfae.inference()` | NMF Autoencoder + inference |
| `nmfkc.rank()` | Rank selection via elbow, cross-validation, ECV, and CPCC |
| `nmfkc.inference()` | Sandwich SE and wild bootstrap p-values for `nmfkc` |
| `nmfkc.DOT()` / `nmfkc.ar.DOT()` / `nmf.sem.DOT()` / `nmfae.DOT()` | Graphviz path diagrams; render with `plot()` |

S3 methods `coef()`, `fitted()`, `residuals()`, `plot()`, `summary()`, `predict()` are available for all model classes. See `?nmfkc` or `browseVignettes("nmfkc")` for the full function list.

# References

- Satoh, K. (2024). Applying Non-negative Matrix Factorization with Covariates to the Longitudinal Data as Growth Curve Model. arXiv:2403.05359. <https://arxiv.org/abs/2403.05359>
- Satoh, K. (2025). Applying non-negative Matrix Factorization with Covariates to Multivariate Time Series Data as a Vector Autoregression Model. *Japanese Journal of Statistics and Data Science*. <https://doi.org/10.1007/s42081-025-00314-0>
- Satoh, K. (2025). Applying non-negative matrix factorization with covariates to label matrix for classification. arXiv:2510.10375. <https://arxiv.org/abs/2510.10375>
- Satoh, K. (2025). Applying non-negative matrix factorization with covariates to structural equation modeling for blind input-output analysis. arXiv:2512.18250. <https://arxiv.org/abs/2512.18250>
- Satoh, K. (2026). Wild Bootstrap Inference for Non-Negative Matrix Factorization with Random Effects. arXiv:2603.01468. <https://arxiv.org/abs/2603.01468>
