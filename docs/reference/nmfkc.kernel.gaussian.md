# Create a Gaussian kernel matrix from covariates

`nmfkc.kernel.gaussian` constructs a Gaussian (RBF) kernel matrix from
covariate matrices. The kernel is defined as \\K(u,v) = \exp(-\beta
\\u - v\\^2)\\. When `V` contains `NA` values, two methods are available
via `na.method`:

- `"pds"`:

  Partial Distance Strategy. Computes the kernel using only observed
  (non-NA) rows, with beta adjusted by \\\beta\_{adj} = \beta \times K /
  K\_{obs}\\ where \\K\\ is the total number of rows and \\K\_{obs}\\ is
  the number of observed rows.

- `"egk"`:

  Expected Gaussian Kernel (Mesquita et al., 2019). Uses a Gaussian
  Mixture Model (GMM) to estimate the conditional distribution of
  missing values given observed values, then computes the expected
  kernel value via a Gamma approximation. Requires `gmm.means`,
  `gmm.sigmas`, and `gmm.weights` passed through `...`.

## Usage

``` r
nmfkc.kernel.gaussian(
  U,
  V = NULL,
  beta = 0.5,
  na.method = c("pds", "egk"),
  ...
)
```

## Source

Mesquita, D., Gomes, J. P., & Rodrigues, L. R. (2019). Gaussian kernels
for incomplete data. *Applied Soft Computing*, 77, 356–365.

## Arguments

- U:

  Covariate matrix \\U(K,N) = (u_1, \dots, u_N)\\. Each row may be
  normalized in advance.

- V:

  Covariate matrix \\V(K,M) = (v_1, \dots, v_M)\\, typically used for
  prediction. If `NULL`, the default is `U`. May contain `NA` values.

- beta:

  Bandwidth parameter for the Gaussian kernel. Default is `0.5`.

- na.method:

  Method for handling `NA` values in `V`. Either `"pds"` or `"egk"`.
  Ignored if `V` has no `NA`.

- ...:

  Additional arguments for EGK method:

  `gmm.G`

  :   Number of GMM components for EGK. Default is `3` (Mesquita et al.,
      2019).

## Value

Kernel matrix \\A(N,M)\\.

## See also

[`nmfkc.kernel`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.md),
[`nmfkc.kernel.beta.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.beta.cv.md),
[`nmfkc.kernel.beta.nearest.med`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.beta.nearest.med.md)

## Examples

``` r
U <- matrix(c(5,10,15,20,25),nrow=1)
V <- matrix(1:25,nrow=1)
A <- nmfkc.kernel.gaussian(U,V,beta=28/1000)
dim(A)
#> [1]  5 25

# PDS example: V with NA in first row
U2 <- matrix(rnorm(20), nrow=2)
V2 <- matrix(rnorm(10), nrow=2)
V2[1, c(2,4)] <- NA
A2 <- nmfkc.kernel.gaussian(U2, V2, beta=0.5, na.method="pds")
```
