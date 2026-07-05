# Methods With Complex Linkage

This documentation covers a few complex methods. In particular:

- L-PLS - Partial Least Squares in L configuration
  ([`lpls`](https://khliland.github.io/multiblock/reference/lpls.md))

- SO-PLS-PM - Sequential and Orthogonalised PLS Path Modeling
  ([`sopls_pm`](https://khliland.github.io/multiblock/reference/SO_TDI.md))

## See also

Overviews of available methods,
[`multiblock`](https://khliland.github.io/multiblock/reference/multiblock.md),
and methods organised by main structure:
[`basic`](https://khliland.github.io/multiblock/reference/basic.md),
[`unsupervised`](https://khliland.github.io/multiblock/reference/unsupervised.md),
[`asca`](https://khliland.github.io/HDANOVA/reference/asca.html),
[`supervised`](https://khliland.github.io/multiblock/reference/supervised.md)
and `complex`.

## Examples

``` r
# L-PLS
sim <- lplsData(I = 30, N = 20, J = 5, K = 6, ncomp = 2)
X1  <- sim$X1; X2 <- sim$X2; X3 <- sim$X3
lp  <- lpls(X1,X2,X3) # exo-L-PLS
```
