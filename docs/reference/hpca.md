# Hierarchical Principal component analysis - HPCA

This is a wrapper for the
[`RGCCA::rgcca`](https://rgcca-factory.github.io/RGCCA/reference/rgcca.html)
function for computing HPCA.

## Usage

``` r
hpca(X, ncomp = 2, scale = FALSE, verbose = FALSE, ...)
```

## Arguments

- X:

  `list` of input blocks.

- ncomp:

  `integer` number of components to extract.

- scale:

  `logical` indicating if variables should be scaled.

- verbose:

  `logical` indicating if diagnostic information should be printed.

- ...:

  additional arguments for RGCCA.

## Value

`multiblock` object including relevant scores and loadings. Relevant
plotting functions:
[`multiblock_plots`](https://khliland.github.io/multiblock/reference/multiblock_plots.md)
and result functions:
[`multiblock_results`](https://khliland.github.io/multiblock/reference/multiblock_results.md).

## Details

HPCA is a hierarchical PCA analysis which combines two or more blocks
into a two-level decomposition with block-wise loadings and scores and
superlevel common loadings and scores. The method is closely related to
the supervised method MB-PLS in structure.

## References

Westerhuis, J.A., Kourti, T., and MacGregor,J.F. (1998). Analysis of
multiblock and hierarchical PCA and PLS models. Journal of Chemometrics,
12, 301–321.

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
pot.hpca   <- hpca(potList)
plot(scores(pot.hpca), labels="names")

```
