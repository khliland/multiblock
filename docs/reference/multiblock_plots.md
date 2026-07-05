# Plot Functions for Multiblock Objects

Plotting procedures for `multiblock` objects.

## Usage

``` r
# S3 method for class 'multiblock'
scoreplot(
  object,
  comps = 1:2,
  block = 0,
  labels,
  identify = FALSE,
  type = "p",
  xlab,
  ylab,
  main,
  ...
)

# S3 method for class 'multiblock'
loadingplot(
  object,
  comps = 1:2,
  block = 0,
  scatter = TRUE,
  labels,
  identify = FALSE,
  type,
  lty,
  lwd = NULL,
  pch,
  cex = NULL,
  col,
  legendpos,
  xlab,
  ylab,
  main,
  pretty.xlabels = TRUE,
  xlim,
  ...
)

loadingweightplot(object, main = "Loading weights", ...)

# S3 method for class 'multiblock'
biplot(
  x,
  block = 0,
  comps = 1:2,
  which = c("x", "y", "scores", "loadings"),
  var.axes = FALSE,
  xlabs,
  ylabs,
  main,
  ...
)

corrplot(object, ...)

# Default S3 method
corrplot(object, ...)

# S3 method for class 'mvr'
corrplot(object, ...)

# S3 method for class 'multiblock'
corrplot(
  object,
  comps = 1:2,
  labels = TRUE,
  col = 1:5,
  plotx = TRUE,
  ploty = TRUE,
  blockScores = FALSE,
  ...
)
```

## Arguments

- object:

  `multiblock` object.

- comps:

  `integer` vector giving components, within block, to plot.

- block:

  `integer/character` for block selection.

- labels:

  `character` indicating if "names" or "numbers" should be plot symbols
  (optional).

- identify:

  `logical` for activating `identify` to interactively identify points.

- type:

  `character` for selecting type of plot to make. Defaults to "p"
  (points) for scatter plots and "l" (lines) for line plots.

- xlab:

  `character` text for x labels.

- ylab:

  `character` text for y labels.

- main:

  `character` text for main header.

- ...:

  Not implemented.

- scatter:

  `logical` indicating if a scatterplot of loadings should be made
  (default = TRUE).

- lty:

  Vector of line type specifications (see
  [`par`](https://rdrr.io/r/graphics/par.html) for details).

- lwd:

  `numeric` vector of line width specifications.

- pch:

  Vector of point specifications (see
  [`points`](https://rdrr.io/r/graphics/points.html) for details).

- cex:

  `numeric` vector of plot size expansions (see
  [`par`](https://rdrr.io/r/graphics/par.html) for details).

- col:

  `integer` vector of symbol/line colours (see
  [`par`](https://rdrr.io/r/graphics/par.html) for details).

- legendpos:

  `character` indicating legend position (if `scatter` is FALSE), e.g.
  `legendpos = "topright"`.

- pretty.xlabels:

  `logical` indicating if xlabels should be more nicely plotted (default
  = TRUE).

- xlim:

  `numeric` vector of length two, with the x limits of the plot
  (optional).

- x:

  `multiblock` object.

- which:

  `character` for selecting type of biplot ("x" = default, "y",
  "scores", "loadings").

- var.axes:

  `logical` indicating if second axes of a biplot should have arrows.

- xlabs:

  `character` vector for labelling first set of biplot points
  (optional).

- ylabs:

  `character` vector for labelling second set of biplot points
  (optional).

- plotx:

  `locical` or `integer`/`character`. Whether to plot the \\X\\
  correlation loadings, optionally which block(s). Defaults to `TRUE`.

- ploty:

  `logical`. Whether to plot the \\Y\\ correlation loadings. Defaults to
  `TRUE`.

- blockScores:

  `logical`. Correlation loadings from blockScores (default = FALSE).

## Value

These plotting routines only generate plots and return no values.

## Details

Plot functions for `scores`, `loadings` and `loading.weights` based on
the functions found in the `pls` package.

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
Common functions for computation and extraction of results are found in
[`multiblock_results`](https://khliland.github.io/multiblock/reference/multiblock_results.md).

## Examples

``` r
data(wine)
sc <- sca(wine[c('Smell at rest', 'View', 'Smell after shaking')], ncomp = 4)
loadingplot(sc, block = 1, labels = "names", scatter = TRUE)

scoreplot(sc, labels = "names")

corrplot(sc)


data(potato)
so <- sopls(Sensory ~ NIRraw + Chemical + Compression, data=potato, ncomp = c(2,2,2), 
            max_comps = 6, validation = "CV", segments = 10)
scoreplot(so, ncomp = c(2,1), block = 3, labels = "names")

corrplot(pcp(so, ncomp = c(2,2,2)))

```
