# Print method for nmfkc.DOT objects

Writes the Graphviz source out as text. Without it, printing a DOT
object fell to `print.default`, which showed the whole graph as one
escaped string (`[1] "digraph NMF \{\n ..."`) — unreadable. Shared by
every `*.DOT` function through the common `"nmfkc.DOT"` class.

## Usage

``` r
# S3 method for class 'nmfkc.DOT'
print(x, ...)
```

## Arguments

- x:

  An object of class `"nmfkc.DOT"` (or a subclass).

- ...:

  Not used.

## Value

`x`, invisibly.

## See also

[`plot.nmfkc.DOT`](https://ksatohds.github.io/nmfkc/reference/plot.nmfkc.DOT.md),
[`nmfkc.DOT`](https://ksatohds.github.io/nmfkc/reference/nmfkc.DOT.md)
