# Canonical Correlation Analysis - CCA

This is a wrapper for the
[`stats::cancor`](https://rdrr.io/r/stats/cancor.html) function for
computing CCA.

## Usage

``` r
cca(X)
```

## Arguments

- X:

  `list` of input data blocks.

## Value

`multiblock` object with associated with printing, scores, loadings.
Relevant plotting functions:
[`multiblock_plots`](https://khliland.github.io/multiblock/reference/multiblock_plots.md)
and result functions:
[`multiblock_results`](https://khliland.github.io/multiblock/reference/multiblock_results.md).

## Details

CCA is a method which maximises correlation between linear combinations
of the columns of two blocks, i.e. max(cor(X1 x a, X2 x b)). This is
done sequentially with deflation in between, such that a sequence of
correlations and weight vectors a and b are associated with a pair of
matrices.

## References

Hotelling, H. (1936) Relations between two sets of variates. Biometrika,
28, 321–377.

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

cca.pot  <- cca(potato[1:2])
```
