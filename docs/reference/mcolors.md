# Colour palette generation from matrix of RGB values

Colour palette generation from matrix of RGB values

## Usage

``` r
mcolors(
  n,
  colmatrix = matrix(c(0, 0, 1, 1, 1, 1, 1, 0, 0), 3, 3, byrow = TRUE)
)
```

## Arguments

- n:

  Integer number of colorus to produce.

- colmatrix:

  A numeric `matrix` of three columns (R,G,B) to generate colour palette
  from.

## Value

A vector of n colours in hexadecimal RGB format.

## Examples

``` r
mcolors(5)
#> [1] "#0000FF" "#8080FF" "#FFFFFF" "#FF8080" "#FF0000"
```
