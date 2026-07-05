# Generalized Canonical Analysis - GCA

This is an interface to both SVD-based (default) and RGCCA-based GCA
(wrapping the `RGCCA::rgcca` function)

## Usage

``` r
gca(X, ncomp = "max", svd = TRUE, tol = 10^-12, corrs = TRUE, ...)
```

## Arguments

- X:

  `list` of input blocks.

- ncomp:

  `integer` number of components to extract, either single integer
  (equal for all blocks), vector (individual per block) or 'max' for
  maximum possible number of components.

- svd:

  `logical` indicating if Singular Value Decomposition approach should
  be used (default=TRUE).

- tol:

  `numeric` tolerance for component inclusion (singular values).

- corrs:

  `logical` indicating if correlations should be calculated for RGCCA
  based approach.

- ...:

  additional arguments for RGCCA approach.

## Value

`multiblock` object including relevant scores and loadings. Relevant
plotting functions:
[`multiblock_plots`](https://khliland.github.io/multiblock/reference/multiblock_plots.md)
and result functions:
[`multiblock_results`](https://khliland.github.io/multiblock/reference/multiblock_results.md).
`blockCoef` contains canonical coefficients, while `blockDecomp`
contains decompositions of each block.

## Details

GCA is a generalisation of Canonical Correlation Analysis to handle
three or more blocks. There are several ways to generalise, and two of
these are available through `gca`. The default is an SVD based approach
estimating a common subspace and measuring mean squared correlation to
this. An alternative approach is available through RGCCA. For the SVD
based approach, the `ncomp` parameter controls the block-wise
decomposition while the following the consensus decomposition is limited
to the minimum number of components from the individual blocks.

## References

- Carroll, J. D. (1968). Generalization of canonical correlation
  analysis to three or more sets of variables. Proceedings of the
  American Psychological Association, pages 227-22.

- Van der Burg, E. and Dijksterhuis, G. (1996). Generalised canonical
  analysis of individual sensory profiles and instrument data, Elsevier,
  pp. 221–258.

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
pot.gca <- gca(potList)
plot(scores(pot.gca), labels="names")

```
