# Structuration des Tableaux à Trois Indices de la Statistique - STATIS

This is a wrapper for the
[`ade4::statis`](https://adeverse.github.io/ade4/reference/statis.html)
function for computing STATIS.

## Usage

``` r
statis(X, ncomp = 3, scannf = FALSE, tol = 1e-07, ...)
```

## Arguments

- X:

  `list` of input blocks.

- ncomp:

  `integer` number of components to extract.

- scannf:

  `logical` indicating if eigenvalue bar plot shoulde be displayed.

- tol:

  `numeric` eigenvalue threshold tolerance.

- ...:

  additional arguments (not used).

## Value

`multiblock` object including relevant scores and loadings. Relevant
plotting functions:
[`multiblock_plots`](https://khliland.github.io/multiblock/reference/multiblock_plots.md)
and result functions:
[`multiblock_results`](https://khliland.github.io/multiblock/reference/multiblock_results.md).

## Details

STATIS is a method, related to MFA, for analysing two or more blocks. It
also decomposes the data into a low-dimensional subspace but uses a
different scaling of the individual blocks.

## References

Lavit, C.; Escoufier, Y.; Sabatier, R.; Traissac, P. (1994). The ACT
(STATIS method). Computational Statistics & Data Analysis. 18: 97

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
are found in
[`multiblock_results`](https://khliland.github.io/multiblock/reference/multiblock_results.md)
and
[`multiblock_plots`](https://khliland.github.io/multiblock/reference/multiblock_plots.md),
respectively.

## Examples

``` r
data(candies)
candyList <- lapply(1:nlevels(candies$candy),function(x)candies$assessment[candies$candy==x,])
can.statis <- statis(candyList)
plot(scores(can.statis), labels="names")

```
