# Generalized Procrustes Analysis - GPA

This is a wrapper for the
[`FactoMineR::GPA`](https://rdrr.io/pkg/FactoMineR/man/gpa.html)
function for computing GPA.

## Usage

``` r
gpa(X, graph = FALSE, ...)
```

## Arguments

- X:

  `list` of input blocks.

- graph:

  `logical` indicating if decomposition should be plotted.

- ...:

  additional arguments for RGCCA approach.

## Value

`multiblock` object including relevant scores and loadings. Relevant
plotting functions:
[`multiblock_plots`](https://khliland.github.io/multiblock/reference/multiblock_plots.md)
and result functions:
[`multiblock_results`](https://khliland.github.io/multiblock/reference/multiblock_results.md).

## Details

GPA is a generalisation of Procrustes analysis, where one matrix is
scaled and rotated to be as similar as possible to another one. Through
the generalisation, individual scaling and rotation of each input matrix
is performed against a common representation which is estimated in an
iterative manner.

## References

Gower, J. C. (1975). Generalized procrustes analysis. Psychometrika. 40:
33–51.

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
data(potato)
potList <- as.list(potato[c(1,2,9)])
pot.gpa    <- gpa(potList)
plot(scores(pot.gpa), labels="names")

```
