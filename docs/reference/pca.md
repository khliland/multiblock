# Principal Component Analysis - PCA

This is a wrapper for the `pls::PCR` function for computing PCA.

## Usage

``` r
pca(X, scale = FALSE, ncomp = 1, ...)
```

## Arguments

- X:

  `matrix` of input data.

- scale:

  `logical` indicating if variables should be standardised
  (default=FALSE).

- ncomp:

  `integer` number of principal components to return.

- ...:

  additional arguments to `pls:pcr`.

## Value

`multiblock` object with scores, loadings, mean X values and explained
variances. Relevant plotting functions:
[`multiblock_plots`](https://khliland.github.io/multiblock/reference/multiblock_plots.md)
and result functions:
[`multiblock_results`](https://khliland.github.io/multiblock/reference/multiblock_results.md).

## Details

PCA is a method for decomposing a matrix into subspace components with
sample scores and variable loadings. It can be formulated in various
ways, but the standard formulation uses singular value decomposition to
create scores and loadings. PCA is guaranteed to be the optimal way of
extracting orthogonal subspaces from a matrix with regard to the amount
of explained variance per component.

## References

Pearson, K. (1901) On lines and planes of closest fit to points in
space. Philosophical Magazine, 2, 559–572.

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
X <- potato$Chemical

pca.pot  <- pca(X, ncomp = 2)
```
