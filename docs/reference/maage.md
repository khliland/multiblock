# Måge plot

Måge plot for SO-PLS
([`sopls`](https://khliland.github.io/multiblock/reference/sopls.md))
cross-validation visualisation.

## Usage

``` r
maage(
  object,
  expl_var = TRUE,
  pure.trace = FALSE,
  pch = 20,
  xlab = "# components",
  ylab = ifelse(expl_var, "Explained variance (%)", "RMSECV"),
  xlim = NULL,
  ylim = NULL,
  cex.text = 0.8,
  ...
)

maageSeq(
  object,
  compSeq = TRUE,
  expl_var = TRUE,
  pch = 20,
  xlab = "# components",
  ylab = ifelse(expl_var, "Explained variance (%)", "RMSECV"),
  xlim = NULL,
  ylim = NULL,
  cex.text = 0.8,
  col = "gray",
  col.block = c("red", "blue", "darkgreen", "purple", "black", "red", "blue",
    "darkgreen"),
  ...
)
```

## Arguments

- object:

  An SO-PLS model (`sopls` object)

- expl_var:

  Logical indicating if explained variance (default) or RMSECV should be
  displayed.

- pure.trace:

  Logical indicating if single block solutions should be traced in the
  plot.

- pch:

  Scalar or symbol giving plot symbol.

- xlab:

  Label for x-axis.

- ylab:

  Label for y-axis.

- xlim:

  Plot limits for x-axis (numeric vector).

- ylim:

  Plot limits for y-axis (numeric vector).

- cex.text:

  Text scaling (scalar) for better readability of plots.

- ...:

  Additional arguments to `plot`.

- compSeq:

  Integer vector giving the sequence of previous components chosen for
  `maageSeq` (see example).

- col:

  Line colour in plot.

- col.block:

  Line colours for blocks (default =
  c('red','blue','darkgreen','purple','black'))

## Value

The `maage` plot has no return.

## Details

This function can either be used for global optimisation across blocks
or sequential optimisation, using `maageSeq`. The examples below show
typical usage.

## See also

Overviews of available methods,
[`multiblock`](https://khliland.github.io/multiblock/reference/multiblock.md),
and methods organised by main structure:
[`basic`](https://khliland.github.io/multiblock/reference/basic.md),
[`unsupervised`](https://khliland.github.io/multiblock/reference/unsupervised.md),
[`asca`](https://khliland.github.io/HDANOVA/reference/asca.html),
[`supervised`](https://khliland.github.io/multiblock/reference/supervised.md)
and
[`complex`](https://khliland.github.io/multiblock/reference/complex.md).

## Examples

``` r
data(wine)
ncomp <- unlist(lapply(wine, ncol))[-5]
so.wine <- sopls(`Global quality` ~ ., data=wine, ncomp=ncomp, 
            max_comps=10, validation="CV", segments=10)
maage(so.wine)


# Sequential search for optimal number of components per block
old.par <- par(mfrow=c(2,2), mar=c(3,3,0.5,1), mgp=c(2,0.7,0))
maageSeq(so.wine)
maageSeq(so.wine, 2)
maageSeq(so.wine, c(2,1))
maageSeq(so.wine, c(2,1,1))

par(old.par)
```
