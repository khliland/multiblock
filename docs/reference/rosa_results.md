# Result functions for ROSA models

Standard result computation and extraction functions for ROSA
([`rosa`](https://khliland.github.io/multiblock/reference/rosa.md)).

## Usage

``` r
# S3 method for class 'rosa'
predict(
  object,
  newdata,
  ncomp = 1:object$ncomp,
  comps,
  type = c("response", "scores"),
  na.action = na.pass,
  ...
)

# S3 method for class 'rosa'
coef(object, ncomp = object$ncomp, comps, intercept = FALSE, ...)

# S3 method for class 'rosa'
print(x, ...)

# S3 method for class 'rosa'
summary(
  object,
  what = c("all", "validation", "training"),
  digits = 4,
  print.gap = 2,
  ...
)

blockexpl(object, ncomp = object$ncomp, type = c("train", "CV"))

# S3 method for class 'rosaexpl'
print(x, digits = 3, compwise = FALSE, ...)

rosa.classify(object, classes, newdata, ncomp, LQ)

# S3 method for class 'rosa'
scores(object, ...)

# S3 method for class 'rosa'
loadings(object, ...)
```

## Arguments

- object:

  A `rosa` object.

- newdata:

  Optional new data with the same types of predictor blocks as the ones
  used for fitting the object.

- ncomp:

  An `integer` giving the number of components to apply (cummulative).

- comps:

  An `integer` vector giving the exact components to apply (subset).

- type:

  For `blockexpl`: Character indicating which type of explained variance
  to compute (default = "train", alternative = "CV").

- na.action:

  Function determining what to do with missing values in `newdata`.

- ...:

  Additional arguments. Currently not implemented.

- intercept:

  A `logical` indicating if coefficients for the intercept should be
  included (default = FALSE).

- x:

  A `rosa` object.

- what:

  A `character` indicating if summary should include all, validation or
  training.

- digits:

  The number of digits used for printing.

- print.gap:

  Gap between columns when printing.

- compwise:

  Logical indicating if block-wise (default/FALSE) or component-wise
  (TRUE) explained variance should be printed.

- classes:

  A `character` vector of class labels.

- LQ:

  A `character` indicating if 'max' (maximum score value), 'lda' or
  'qda' should be used when classifying.

## Value

Returns depend on method used, e.g. `predict.rosa` returns predicted
responses or scores depending on inputs, `coef.rosa` returns regression
coefficients, `blockexpl` returns an object of class `rosaexpl`
containing block-wise and component-wise explained variance contained in
a matrix with attributes.

## Details

Usage of the functions are shown using generics in the examples below.
Prediction, regression coefficients, object printing and summary are
available through: `predict.rosa`, `coef.rosa`, `print.rosa` and
`summary.rosa`. Explained variances are available (block-wise and
global) through `blockexpl` and `print.rosaexpl`. Scores and loadings
have their own extensions of
[`scores()`](https://khliland.github.io/pls/reference/scores.html) and
[`loadings()`](https://khliland.github.io/pls/reference/scores.html)
throught `scores.rosa` and `loadings.rosa`. Finally, there is work in
progress on classifcation support through `rosa.classify`.

When `type` is `"response"` (default), predicted response values are
returned. If `comps` is missing (or is `NULL`), predictions for
`length(ncomp)` models with `ncomp[1]` components, `ncomp[2]`
components, etc., are returned. Otherwise, predictions for a single
model with the exact components in `comps` are returned. (Note that in
both cases, the intercept is always included in the predictions. It can
be removed by subtracting the `Ymeans` component of the fitted model.)

If `comps` is missing (or is `NULL`), `coef()[,,ncomp[i]]` are the
coefficients for models with `ncomp[i]` components, for \\i = 1, \ldots,
length(ncomp)\\. Also, if `intercept = TRUE`, the first dimension is
\\nxvar + 1\\, with the intercept coefficients as the first row.

If `comps` is given, however, `coef()[,,comps[i]]` are the coefficients
for a model with only the component `comps[i]`, i.e., the contribution
of the component `comps[i]` on the regression coefficients.

## References

Liland, K.H., Næs, T., and Indahl, U.G. (2016). ROSA - a fast extension
of partial least squares regression for multiblock data analysis.
Journal of Chemometrics, 30, 651–662, doi:10.1002/cem.2824.

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
Common functions for computation and extraction of results and plotting
are found in `rosa_results` and
[`rosa_plots`](https://khliland.github.io/multiblock/reference/rosa_plots.md),
respectively.

## Examples

``` r
data(potato)
mod <- rosa(Sensory[,1] ~ ., data = potato, ncomp = 5, subset = 1:20)
testset <- potato[-(1:20),]; testset$Sensory <- testset$Sensory[,1,drop=FALSE]
predict(mod, testset, ncomp=5)
#> , , 5 comps
#> 
#>    Sensory[, 1]
#> 21     4.892899
#> 22     4.454407
#> 23     8.177141
#> 24     3.087408
#> 25     4.671565
#> 26     3.926871
#> 
dim(coef(mod, ncomp=5)) # <variables x responses x components>
#> [1] 3946    1    1
print(mod)
#> Response Orinented Sequential Alternation , fitted with the CPPLS algorithm.
#> Call:
#> rosa(formula = Sensory[, 1] ~ ., ncomp = 5, data = potato, subset = 1:20)
summary(mod)
#> Data:    X dimension: 20 3946 
#>  Y dimension: 20 1
#> Fit method:
#> Number of components considered: 5
#> TRAINING: % variance explained
#>               1 comps  2 comps  3 comps  4 comps  5 comps
#> X               28.09    45.00    49.16    62.24    72.58
#> Sensory[, 1]    72.84    91.32    93.68    95.15    96.73
blockexpl(mod)
#> Block-wise explained variance
#> 
#>   Chemical Compression NIRraw NIRcooked CPMGraw CPMGcooked FIDraw FIDcooked
#> X    0.453       0.103  0.169         0       0          0      0         0
#> Y    0.767       0.016  0.185         0       0          0      0         0
#>   residual
#> X    0.274
#> Y    0.033
print(blockexpl(mod), compwise=TRUE)
#> Component-wise explained variance
#> 
#>   comp.1 (Chemical) comp.2 (NIRraw) comp.3 (Chemical) comp.4 (Chemical)
#> X             0.281           0.169             0.042             0.131
#> Y             0.728           0.185             0.024             0.015
#>   comp.5 (Compression) residual
#> X                0.103    0.033
#> Y                0.016    0.033
```
