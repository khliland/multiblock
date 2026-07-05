# L-PLS data simulation for exo-type analysis

Three data blocks are simulated to express covariance in an exo-L-PLS
direction (see
[`lpls`](https://khliland.github.io/multiblock/reference/lpls.md).
Dimensionality and number of underlying components can be controlled.

## Usage

``` r
lplsData(I = 30, N = 20, J = 5, K = 6, ncomp = 2)
```

## Arguments

- I:

  `numeric` number of rows of X1 and X2

- N:

  `numeric` number of columns in X1 and X3

- J:

  `numeric` number of columns in X2

- K:

  `numeric` number of rows in X3

- ncomp:

  `numeric` number of latent components

## Value

A `list` of three matrices with dimensions matching in an L-shape.

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

## Author

Solve Sæbø (adapted by Kristian Hovde Liland)

## Examples

``` r
lp <- lplsData(I = 30, N = 20, J = 5, K = 6, ncomp = 2)
names(lp)
#> [1] "X1" "X2" "X3"
```
