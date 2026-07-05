# Result functions for L-PLS objects (`lpls`)

Correlation loading plot, prediction and cross-validation for L-PLS
models with class
[`lpls`](https://khliland.github.io/multiblock/reference/lpls.md).

## Usage

``` r
# S3 method for class 'lpls'
plot(
  x,
  comps = c(1, 2),
  doplot = c(TRUE, TRUE, TRUE),
  level = c(2, 2, 2),
  arrow = c(1, 0, 1),
  xlim = c(-1, 1),
  ylim = c(-1, 1),
  samplecol = 4,
  pathcol = 2,
  varcol = "grey70",
  varsize = 1,
  sampleindex = 1:dim(x$corloadings$R22)[1],
  pathindex = 1:dim(x$corloadings$R3)[1],
  varindex = 1:dim(x$corloadings$R21)[1],
  ...
)

# S3 method for class 'lpls'
predict(
  object,
  X1new = NULL,
  X2new = NULL,
  X3new = NULL,
  exo.direction = c("X2", "X3"),
  ...
)

lplsCV(object, segments1 = NULL, segments2 = NULL, trace = TRUE)
```

## Arguments

- x:

  `lpls` object

- comps:

  `integer` vector of components.

- doplot:

  `logical` indicating if plotting should be performed.

- level:

  `integer` vector of length 3 for selecting plot symbol. 1=dots.
  2=dimnames.

- arrow:

  `integer` vector of length 3 indicating arrows (1) or not (0).

- xlim:

  `numeric` x limits.

- ylim:

  `numeric` y limits.

- samplecol:

  `character` for sample colours.

- pathcol:

  `character` for third colour.

- varcol:

  `character` for variable colours.

- varsize:

  `numeric` size of symbols for variables.

- sampleindex:

  `integer` for selecting samples.

- pathindex:

  `integer` for selecting in third direction.

- varindex:

  `integer` for selecting variables.

- ...:

  Not implemented.

- object:

  `lpls` object.

- X1new:

  `matrix` of new X1 samples.

- X2new:

  `matrix` of new X2 samples.

- X3new:

  `matrix` of new X3 samples.

- exo.direction:

  `character` selecting "X2" or "X3" prediction.

- segments1:

  `list` of sample segments.

- segments2:

  `list` of variable segments.

- trace:

  `logical` indicating if verbose mode should be selected.

## Value

Nothing is return for plotting (`plot.lpls`), predicted values are
returned for predictions (`predict.lpls`) and cross-validation metrics
are returned for for cross-validation (`lplsCV`).

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
# Simulate data set
sim <- lplsData(I = 30, N = 20, J = 5, K = 6, ncomp = 2)
X1  <- sim$X1; X2 <- sim$X2; X3 <- sim$X3

# exo-L-PLS:
lp.exo  <- lpls(X1,X2,X3, ncomp = 2)
# Predict X1
pred.exo.X2 <- predict(lp.exo, X1new = X1, exo.direction = "X2")
# Predict X3
pred.exo.X2 <- predict(lp.exo, X1new = X1, exo.direction = "X3")

# endo-L-PLS:
lp.endo <- lpls(X1,X2,X3, ncomp = 2, type = "endo")
# Predict X1 from X2 and X3 (in this case fitted values):
pred.endo.X1 <- predict(lp.endo, X2new = X2, X3new = X3)

# LOO cross-validation horizontally
lp.cv1 <- lplsCV(lp.exo, segments1 = as.list(1:dim(X1)[1]))
#> Segment 1 of 30 completed
#> Segment 2 of 30 completed
#> Segment 3 of 30 completed
#> Segment 4 of 30 completed
#> Segment 5 of 30 completed
#> Segment 6 of 30 completed
#> Segment 7 of 30 completed
#> Segment 8 of 30 completed
#> Segment 9 of 30 completed
#> Segment 10 of 30 completed
#> Segment 11 of 30 completed
#> Segment 12 of 30 completed
#> Segment 13 of 30 completed
#> Segment 14 of 30 completed
#> Segment 15 of 30 completed
#> Segment 16 of 30 completed
#> Segment 17 of 30 completed
#> Segment 18 of 30 completed
#> Segment 19 of 30 completed
#> Segment 20 of 30 completed
#> Segment 21 of 30 completed
#> Segment 22 of 30 completed
#> Segment 23 of 30 completed
#> Segment 24 of 30 completed
#> Segment 25 of 30 completed
#> Segment 26 of 30 completed
#> Segment 27 of 30 completed
#> Segment 28 of 30 completed
#> Segment 29 of 30 completed
#> Segment 30 of 30 completed

# LOO cross-validation vertically
lp.cv2 <- lplsCV(lp.exo, segments2 = as.list(1:dim(X1)[2]))
#> Segment 1 of 20 completed
#> Segment 2 of 20 completed
#> Segment 3 of 20 completed
#> Segment 4 of 20 completed
#> Segment 5 of 20 completed
#> Segment 6 of 20 completed
#> Segment 7 of 20 completed
#> Segment 8 of 20 completed
#> Segment 9 of 20 completed
#> Segment 10 of 20 completed
#> Segment 11 of 20 completed
#> Segment 12 of 20 completed
#> Segment 13 of 20 completed
#> Segment 14 of 20 completed
#> Segment 15 of 20 completed
#> Segment 16 of 20 completed
#> Segment 17 of 20 completed
#> Segment 18 of 20 completed
#> Segment 19 of 20 completed
#> Segment 20 of 20 completed

# Three-fold CV, horizontal
lp.cv3 <- lplsCV(lp.exo, segments1 = as.list(1:10, 11:20, 21:30))
#> Segment 1 of 10 completed
#> Segment 2 of 10 completed
#> Segment 3 of 10 completed
#> Segment 4 of 10 completed
#> Segment 5 of 10 completed
#> Segment 6 of 10 completed
#> Segment 7 of 10 completed
#> Segment 8 of 10 completed
#> Segment 9 of 10 completed
#> Segment 10 of 10 completed
```
