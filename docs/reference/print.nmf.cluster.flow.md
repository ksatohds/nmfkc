# Print method for nmf.cluster.flow objects

Prints a one-line header, the adjacent-result ARI, and the \\N \times
R\\ cluster table (rows = individuals, columns = results, entries =
cluster number). Use
[`plot`](https://rdrr.io/r/graphics/plot.default.html) for the diagram.

## Usage

``` r
# S3 method for class 'nmf.cluster.flow'
print(x, ...)
```

## Arguments

- x:

  An object of class `"nmf.cluster.flow"`.

- ...:

  Passed to the table's `print`.

## Value

`x`, invisibly.
