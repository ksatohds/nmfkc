# Create a kernel matrix from covariates

`nmfkc.kernel` constructs a kernel matrix from covariate matrices. It
supports Gaussian, Exponential, Periodic, Linear, Normalized Linear, and
Polynomial kernels.

## Usage

``` r
nmfkc.kernel(
  U,
  V = NULL,
  kernel = c("Gaussian", "Exponential", "Periodic", "Linear", "NormalizedLinear",
    "Polynomial"),
  ...
)
```

## Source

Satoh, K. (2026). Applying Non-negative Matrix Factorization with
Covariates to the Longitudinal Data as Growth Curve Model. *American
Journal of Mathematical and Management Sciences*. In press.
arXiv:2403.05359. <https://arxiv.org/abs/2403.05359>

## Arguments

- U:

  Covariate matrix \\U(K,N) = (u_1, \dots, u_N)\\. Each row may be
  normalized in advance.

- V:

  Covariate matrix \\V(K,M) = (v_1, \dots, v_M)\\, typically used for
  prediction. If `NULL`, the default is `U`.

- kernel:

  Kernel function to use. Default is `"Gaussian"`. Options are
  `"Gaussian"`, `"Exponential"`, `"Periodic"`, `"Linear"`,
  `"NormalizedLinear"`, and `"Polynomial"`.

- ...:

  Additional arguments passed to the specific kernel function (e.g.,
  `beta`, `degree`).

## Value

Kernel matrix \\A(N,M)\\.

## See also

[`nmfkc.kernel.gaussian`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.gaussian.md),
[`nmfkc.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.cv.md)

## Examples

``` r
# Example.
Y <- matrix(cars$dist,nrow=1)
U <- matrix(c(5,10,15,20,25),nrow=1)
V <- matrix(cars$speed,nrow=1)
A <- nmfkc.kernel(U,V,beta=28/1000)
dim(A)
#> [1]  5 50
result <- nmfkc(Y,A,rank=1)
#> Y(1,50)~X(1,1)C(1,5)A(5,50)=XB(1,50)...
#> 0sec
plot(as.vector(V),as.vector(Y))
lines(as.vector(V),as.vector(result$XB),col=2,lwd=2)
```
