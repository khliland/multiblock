# ASCA Result Methods

Various plotting procedures for
[`asca`](https://khliland.github.io/multiblock/reference/asca.md)
objects.

## Usage

``` r
# S3 method for class 'asca'
loadingplot(object, factor = 1, comps = 1:2, ...)

# S3 method for class 'asca'
scoreplot(
  object,
  factor = 1,
  comps = 1:2,
  pch.scores = 19,
  pch.projections = 1,
  gr.col = 1:nlevels(object$effects[[factor]]),
  ellipsoids,
  confidence,
  xlim,
  ylim,
  xlab,
  ylab,
  legendpos,
  ...
)
```

## Arguments

- object:

  `asca` object.

- factor:

  `integer/character` for selecting a model factor.

- comps:

  `integer` vector of selected components.

- ...:

  additional arguments to underlying methods.

- pch.scores:

  `integer` plotting symbol.

- pch.projections:

  `integer` plotting symbol.

- gr.col:

  `integer` vector of colours for groups.

- ellipsoids:

  `character` "confidence" or "data" ellipsoids for balanced fixed
  effect models.

- confidence:

  `numeric` vector of ellipsoid confidences, default = c(0.4, 0.68,
  0.95).

- xlim:

  `numeric` x limits.

- ylim:

  `numeric` y limits.

- xlab:

  `character` x label.

- ylab:

  `character` y label.

- legendpos:

  `character` position of legend.

## Value

The plotting routines have no return.

## Details

Usage of the functions are shown using generics in the examples in
[`asca`](https://khliland.github.io/multiblock/reference/asca.md). Plot
routines are available as `scoreplot.asca` and `loadingplot.asca`.

## References

- Smilde, A., Jansen, J., Hoefsloot, H., Lamers,R., Van Der Greef, J.,
  and Timmerman, M.(2005). ANOVA-Simultaneous Component Analysis (ASCA):
  A new tool for analyzing designed metabolomics data. Bioinformatics,
  21(13), 3043–3048.

- Liland, K.H., Smilde, A., Marini, F., and Næs,T. (2018). Confidence
  ellipsoids for ASCA models based on multivariate regression theory.
  Journal of Chemometrics, 32(e2990), 1–13.

- Martin, M. and Govaerts, B. (2020). LiMM-PCA: Combining ASCA+ and
  linear mixed models to analyse high-dimensional designed data. Journal
  of Chemometrics, 34(6), e3232.

## See also

Overviews of available methods,
[`multiblock`](https://khliland.github.io/multiblock/reference/multiblock.md),
and methods organised by main structure:
[`basic`](https://khliland.github.io/multiblock/reference/basic.md),
[`unsupervised`](https://khliland.github.io/multiblock/reference/unsupervised.md),
[`asca`](https://khliland.github.io/multiblock/reference/asca.md),
[`supervised`](https://khliland.github.io/multiblock/reference/supervised.md)
and
[`complex`](https://khliland.github.io/multiblock/reference/complex.md).
Common functions for computation and extraction of results are found in
[`asca_results`](https://khliland.github.io/multiblock/reference/asca_results.md).
