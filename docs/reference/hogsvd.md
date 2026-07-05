# Higher Order Generalized SVD - HOGSVD

This is a simple implementation for computing HOGSVD

## Usage

``` r
hogsvd(X)
```

## Arguments

- X:

  `list` of input blocks.

## Value

`multiblock` object including relevant scores and loadings. Relevant
plotting functions:
[`multiblock_plots`](https://khliland.github.io/multiblock/reference/multiblock_plots.md)
and result functions:
[`multiblock_results`](https://khliland.github.io/multiblock/reference/multiblock_results.md).

## Details

HOGSVD is a generalisation of SVD to two or more blocks. It finds a
common set of loadings across blocks and individual sets of scores per
block.

## References

Ponnapalli, S. P., Saunders, M. A., Van Loan, C. F., & Alter, O. (2011).
A higher-order generalized singular value decomposition for comparison
of global mRNA expression from multiple organisms. PloS one, 6(12),
e28072.

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
can.hogsvd <- hogsvd(candyList)
scoreplot(can.hogsvd, block=1, labels="names")

```
