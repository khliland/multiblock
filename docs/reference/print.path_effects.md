# Print Path Effects

Prints a summary of path effects analysis with optional bootstrap
standard errors.

## Usage

``` r
# S3 method for class 'path_effects'
print(
  x,
  from = NULL,
  to = NULL,
  nsmall = 2,
  digits = 1,
  scaled = FALSE,
  individual = FALSE,
  boot = TRUE,
  ...
)
```

## Arguments

- x:

  An object of class `path_effects`.

- from:

  Optional source block name or index. If specified with `to`, prints
  only that block pair.

- to:

  Optional target block name or index. If specified with `from`, prints
  only that block pair.

- nsmall:

  Number of decimal places to display. Default is 2.

- digits:

  Number of significant digits. Default is 1.

- scaled:

  Logical. If `TRUE`, display scaled contributions. Default is `FALSE`.

- individual:

  Logical. If `TRUE`, display individual-level effects. Default is
  `FALSE`.

- boot:

  Logical. If `TRUE` and bootstrap samples are available, display
  standard errors. Default is `TRUE`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`.

## See also

[`path_effects()`](https://khliland.github.io/multiblock/reference/path_effects.md),
[`plot.path_effects()`](https://khliland.github.io/multiblock/reference/plot.path_effects.md)
