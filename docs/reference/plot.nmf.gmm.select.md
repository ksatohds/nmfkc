# Plot method for nmf.gmm.select objects

Draws BIC and ICL against \\K\\, with the selected \\K\\ marked. Every
other selector in the package
([`nmfkc.rank`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rank.md),
[`nmfkc.ard`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ard.md),
[`nmfkc.consensus`](https://ksatohds.github.io/nmfkc/reference/nmfkc.consensus.md),
[`nmfkc.bicv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.bicv.md))
has a plot method; without this one `plot(nmf.gmm.select(...))` fell
through to `plot.default` and errored with "'x' is a list, but does not
have components 'x' and 'y'".

## Usage

``` r
# S3 method for class 'nmf.gmm.select'
plot(x, ...)
```

## Arguments

- x:

  An object of class `"nmf.gmm.select"`.

- ...:

  Passed to the underlying plot call.

## Value

`x`, invisibly.

## Development status

Part of the **experimental** NMF-GMM family: see
[`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md). The
interface may change in future versions — argument names, defaults and
the contents of the returned object are not yet stable.

## See also

[`nmf.gmm.select`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.select.md),
[`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md)
