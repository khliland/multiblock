# Plot Path Effects

Visualisation of multiblock path effects decomposition.

## Usage

``` r
# S3 method for class 'path_effects'
plot(
  x,
  from = NULL,
  to = NULL,
  scaled = FALSE,
  individual = FALSE,
  spectra = numeric(0),
  mar = c(0.5, 2.8, 1.8, 0.3),
  ...
)
```

## Arguments

- x:

  An object of class `path_effects`.

- from:

  Optional source block name or index. If specified with `to`, plots
  only that block pair.

- to:

  Optional target block name or index. If specified with `from`, plots
  only that block pair.

- scaled:

  Logical. If `TRUE`, display scaled contributions. Default is `FALSE`.

- individual:

  Logical. If `TRUE`, display individual-level effects. Default is
  `FALSE`.

- spectra:

  Optional numeric vector for spectral data visualization. Default is
  empty.

- mar:

  Margins for the plot. Default is `c(0.5,2.8,1.8,0.3)`.

- ...:

  Additional arguments passed to plotting functions.

## Value

Invisibly returns `NULL`.

## See also

[`path_effects()`](https://khliland.github.io/multiblock/reference/path_effects.md),
[`print.path_effects()`](https://khliland.github.io/multiblock/reference/print.path_effects.md)
