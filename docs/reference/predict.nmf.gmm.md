# Class assignments / responsibilities from an NMF-GMM fit

Class assignments / responsibilities from an NMF-GMM fit

## Usage

``` r
# S3 method for class 'nmf.gmm'
predict(object, type = c("class", "responsibility"), ...)
```

## Arguments

- object:

  An object of class `"nmf.gmm"`.

- type:

  `"class"` (default) returns the hard cluster labels;
  `"responsibility"` returns the N x K posterior-probability matrix.

- ...:

  Ignored.

## Value

A vector of labels or the responsibility matrix.

## Development status

Part of the **experimental** NMF-GMM family: see
[`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md). The
interface may change in future versions — argument names, defaults and
the contents of the returned object are not yet stable.

## See also

[`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md),
[`nmf.gmm.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.inference.md),
[`nmf.gmm.select`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.select.md)
