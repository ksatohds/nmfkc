# Plot method for nmfkc.DOT objects

Renders a DOT graph string using
[`DiagrammeR::grViz`](https://rich-iannone.github.io/DiagrammeR/reference/grViz.html).
If the DiagrammeR package is not installed, prints the DOT source to the
console instead.

This method handles all DOT objects produced by the nmfkc package:
[`nmfkc.DOT`](https://ksatohds.github.io/nmfkc/reference/nmfkc.DOT.md),
[`nmf.rrr.DOT`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.DOT.md),
[`nmf.ffb.DOT`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.DOT.md),
and
[`nmfkc.ar.DOT`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.DOT.md).

## Usage

``` r
# S3 method for class 'nmfkc.DOT'
plot(x, ...)
```

## Arguments

- x:

  An object of class `"nmfkc.DOT"` (or a subclass thereof).

- ...:

  Not used.

## Value

Called for its side effect (rendering). Returns `x` invisibly.

## See also

[`nmfkc.DOT`](https://ksatohds.github.io/nmfkc/reference/nmfkc.DOT.md),
[`nmf.rrr.DOT`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.DOT.md),
[`nmf.ffb.DOT`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.DOT.md),
[`nmfkc.ar.DOT`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.DOT.md),
[`print.nmfkc.DOT`](https://ksatohds.github.io/nmfkc/reference/print.nmfkc.DOT.md)
