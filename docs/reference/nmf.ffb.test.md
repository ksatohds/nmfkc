# Test the feed-forward null of an NMF-FFB fit

Calibrated test of the hypothesis that the outcomes need no
indicator-level feedback: the feed-forward model \\Y_1 = X(\Theta_2
Y_2 + U) + \mathcal{E}\\ with a full random-effect covariance against
the feedback model with \\\Theta_1 \ne 0\\. This is the only inferential
step of the procedure (see the workflow below);
[`nmf.ffb.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.inference.md)
is for the uncertainty of the coefficients *after* the null has been
rejected, and
[`nmf.ffb.diagnostics`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.diagnostics.md)
for what the selected feedback looks like.

## Usage

``` r
nmf.ffb.test(
  object,
  Y1,
  Y2,
  B = 1000L,
  calibration = c("procedure", "conditional"),
  seed = 123L,
  ...
)
```

## Arguments

- object:

  A fit from
  [`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md)
  with `method = "fiml"`.

- Y1, Y2:

  The endogenous and exogenous matrices the fit was computed from
  (variables in rows, units in columns).

- B:

  Number of bootstrap replicates (default 1000). The smallest reportable
  \\p\\-value is \\1/(1+B)\\.

- calibration:

  What is re-applied to each null replicate; see above. `"conditional"`
  is selected automatically, with a warning, when the fit used a basis
  or an exclusion restriction supplied by the caller, since there is
  then nothing to re-estimate.

- seed:

  Base seed; replicate `b` uses `seed + b`. The caller's random stream
  is restored on exit.

- ...:

  Passed to the bootstrap: `cores` (or `ncores`) to run the replicates
  in parallel, defaulting to `getOption("mc.cores", 1L)` as elsewhere in
  the package; `fiml.maxit` and `factr` to control the optimizer in the
  replicates; `print.trace = TRUE` for progress.

## Value

An object of class `"nmf.ffb.test"`:

- LR, LR.df:

  Observed `c(full, selected)` likelihood ratios and the number of free
  feedback entries left by the exclusion restriction.

- LR.p.boot:

  Calibrated \\p\\-values of the two statistics.

- LR.null.quantile:

  The 0.95 quantiles of the null distributions.

- LR.boot:

  \\B \times 2\\ matrix of the null replicates.

- prob.select.null:

  Null false-selection rate: the probability that BIC retains at least
  one feedback entry when the feed-forward model is true.

- C1.restriction.change.rate:

  (`"procedure"` only) share of null replicates in which the exclusion
  restriction moved.

- LR.boot.df:

  Free entries per null replicate (varies under `"procedure"`).

- LR.boot.n.ok, LR.boot.n.nonconv:

  Usable replicates, and how many did not meet the optimizer tolerance.

## Details

The likelihood ratio does not have a \\\chi^2\\ null distribution: the
constraint \\\Theta_1 \ge 0\\ puts the null on the boundary of the
parameter space, and the support of \\\Theta_1\\ is chosen from the same
data. It is therefore calibrated by a parametric bootstrap from the
fitted feed-forward null, \$\$Y_1^{\*} = (\hat X(\hat\Theta_2 Y_2 +
U^{\*}) + \mathcal{E}^{\*})\_{+}, \quad U^{\*} \sim N(0, \hat\Phi), \\
\mathcal{E}^{\*} \sim N(0, \hat\Psi),\$\$ with negative draws clipped at
zero, and the reported \\p\\-value is \\(1 + \\\\LR^{\*} \ge LR\\) /
(1 + B\_{ok})\\.

### What is re-applied to each replicate

`calibration = "procedure"` (the default) re-applies the *whole
procedure*: stage 1 is re-run on each \\Y_1^{\*}\\, the basis is
re-estimated, the exclusion restriction is re-derived from that basis,
and stage 2 follows. The \\p\\-value is then the operating
characteristic of the complete exploratory procedure, which is what is
applied to a fresh data set. This is the calibration to use.

`"conditional"` holds \\\hat X\\ and the exclusion restriction at their
fitted values and re-runs stage 2 only. It is valid only if basis and
restriction are independent of the \\Y_1\\ being tested, which they are
not when all three come from the same sample; it is provided for
comparison with that practice, and it is anti-conservative.

Sample splitting – basis and restriction from one half of the units,
test on other – also removes the dependence, and was offered here in
0.9.8. It was withdrawn in 0.9.8: measured against `"procedure"` it has
the same size and lower power, because the test uses \\N/2\\ units, and
it cannot be run at all when the halves are too small for stage 1.

### Reading the result

Report `LR.p.boot` together with `prob.select.null` and, for
`calibration = "procedure"`, `C1.restriction.change.rate`. The last two
say whether the exclusion restriction is determined well enough for the
test to mean anything: if the dominant factor of an outcome moves from
one null replicate to the next, a previously blocked entry becomes free
and a feedback coefficient can absorb loading structure that the
feed-forward model attributes to \\X\\. A large `prob.select.null` with
a large null quantile is the signature of a data set in which the
procedure finds feedback whether or not there is any.

## Workflow

    ecv <- nmf.ffb.ecv(Y1, Y2, rank = 1:5)   # 1. choose Q
    fit <- nmf.ffb(Y1, Y2, rank = ecv$rank)  # 2. estimate; BIC selects the support
    tst <- nmf.ffb.test(fit, Y1, Y2)         # 3. test the feed-forward null  <- only inference
    dgn <- nmf.ffb.diagnostics(fit)          # 4. cycles, spectral radius, identifiability
    inf <- nmf.ffb.inference(fit, Y1, Y2)    # 5. intervals for the retained entries

## See also

[`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md),
[`nmf.ffb.ecv`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.ecv.md),
[`nmf.ffb.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.inference.md),
[`nmf.ffb.diagnostics`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.diagnostics.md)
