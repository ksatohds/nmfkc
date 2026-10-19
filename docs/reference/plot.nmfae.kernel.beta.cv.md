# Plot method for nmfae.kernel.beta.cv objects

Displays the cross-validation objective function across candidate `beta`
values (log scale). The optimal beta is highlighted in red.

## Usage

``` r
# S3 method for class 'nmfae.kernel.beta.cv'
plot(x, ...)
```

## Arguments

- x:

  An object of class `"nmfae.kernel.beta.cv"` returned by
  [`nmf.rrr.kernel.beta.cv`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.kernel.beta.cv.md).

- ...:

  Additional graphical parameters passed to `plot`.

## Value

Invisible `NULL`. Called for its side effect (plot).

## See also

[`nmf.rrr.kernel.beta.cv`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.kernel.beta.cv.md)
