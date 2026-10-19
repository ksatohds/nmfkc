# Extract the covariate-coefficient matrix from an NMF-GMM fit

Extract the covariate-coefficient matrix from an NMF-GMM fit

## Usage

``` r
# S3 method for class 'nmf.gmm'
coef(object, ...)
```

## Arguments

- object:

  An object of class `"nmf.gmm"`.

- ...:

  Ignored.

## Value

The \\Q\times R\\ coefficient matrix \\C\\ (\\=\Theta\\).

## Development status

Part of the **experimental** NMF-GMM family: see
[`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md). The
interface may change in future versions — argument names, defaults and
the contents of the returned object are not yet stable.

## See also

[`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md),
[`nmf.gmm.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.inference.md),
[`nmf.gmm.select`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.select.md)
