# Fitted (responsibility-averaged) reconstruction of an NMF-GMM fit

Fitted (responsibility-averaged) reconstruction of an NMF-GMM fit

## Usage

``` r
# S3 method for class 'nmf.gmm'
fitted(object, ...)
```

## Arguments

- object:

  An object of class `"nmf.gmm"`.

- ...:

  Ignored.

## Value

The \\P\times N\\ fitted matrix \\\hat Y = X(CA + \mu\gamma^\top)\\.

## Development status

Part of the **experimental** NMF-GMM family: see
[`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md). The
interface may change in future versions — argument names, defaults and
the contents of the returned object are not yet stable.

## See also

[`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md),
[`nmf.gmm.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.inference.md),
[`nmf.gmm.select`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.select.md)
