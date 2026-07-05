# Analyse Multivariate Path Effects

Decomposes effects between blocks in a multiblock system into total
effects, unique effects, common contributions, and additional effects.
Supports both SO-PLS (Sequential Orthogonalised PLS) and ordinary least
squares (OLS) approaches with cross-validation or fitted values.

## Usage

``` r
path_effects(
  relations,
  blocks,
  validation,
  segments,
  SO = TRUE,
  fits = FALSE,
  boot = 0,
  ncomp = NULL,
  transitive_closure = FALSE,
  ...
)
```

## Arguments

- relations:

  A matrix defining the path structure. Rows specify relationships
  between blocks, with 2 columns: `c(from_block_index, to_block_index)`.

- blocks:

  A multiblock data frame or list of named matrices containing the block
  data.

- validation:

  Validation method passed to
  [`pls::plsr()`](https://khliland.github.io/pls/reference/mvr.html).
  Common options: `"CV"` for cross-validation, "LOO" for leave-one-out
  validation, or `"none"` for fitted values.

- segments:

  Number of segments for cross-validation. Default is 5.

- SO:

  Logical. If `TRUE` (default), use SO-PLS; if `FALSE`, use OLS.

- fits:

  Logical. If `TRUE`, use fitted values instead of cross-validation.
  Default is `FALSE`.

- boot:

  Number of bootstrap samples for computing standard errors. Default is
  0 (no bootstrapping).

- ncomp:

  Optional component settings per predictor block. Supply either:

  - a vector with length equal to the number of predictor blocks
    (`unique(relations[,1])`), in ascending predictor-block order, or

  - a full-length vector with one entry per block (non-predictors may be
    `0`). Positive integers set upper limits; negative integers force
    exact component counts. All predictor entries must be non-zero and
    have the same sign.

- transitive_closure:

  Logical. If `TRUE`, automatically add transitive closure to the
  relations. Default is `FALSE`.

- ...:

  Additional arguments passed to underlying fitting functions.

## Value

An object of class `path_effects`, which is a matrix with the following
components for each path:

- `Tot`: Total effect

- `Un`: Unique effect

- `Co`: Common contribution

- `Ad`: Additional effect

Attributes include:

- `scaled`: Scaled contribution matrix

- `individual`: Individual-level contributions

- `boot`: Bootstrap replicates (if `boot > 0`)

## See also

[`print.path_effects()`](https://khliland.github.io/multiblock/reference/print.path_effects.md),
[`plot.path_effects()`](https://khliland.github.io/multiblock/reference/plot.path_effects.md)

## Examples
