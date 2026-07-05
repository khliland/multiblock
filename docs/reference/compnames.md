# Vector of component names

Convenience function for creating a vector of component names based on
the dimensions the input object (`matrix` or object having a `score`
function).

## Usage

``` r
compnames(object, comps, explvar = FALSE, ...)
```

## Arguments

- object:

  An object fitted using the multiblock package.

- comps:

  `integer` vector of components.

- explvar:

  `logical` indicating if explained variances should be included.

- ...:

  Unused

## Value

A `character` vector of component names.

## Details

This is a copy of `compnames` from the `pls` package to work with
`multiblock` objects.
