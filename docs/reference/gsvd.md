# Generalised Singular Value Decomposition - GSVD

This is a wrapper for the
[`geigen::gsvd`](https://rdrr.io/pkg/geigen/man/gsvd.html) function for
computing GSVD.

## Usage

``` r
gsvd(X)
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

GSVD is a generalisation of SVD to two variable-linked matrices where
common loadings and block-wise scores are estimated.

## References

Van Loan, C. (1976) Generalizing the singular value decomposition. SIAM
Journal on Numerical Analysis, 13, 76–83.

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

gsvd.pot <- gsvd(lapply(potato[3:4], t))
```
