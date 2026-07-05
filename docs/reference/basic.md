# Single- and Two-Block Methods

This documentation covers a range of single- and two-block methods. In
particular:

- PCA - Principal Component Analysis
  ([`pca`](https://khliland.github.io/multiblock/reference/pca.md))

- PCR - Principal Component Regression
  ([`pcr`](https://khliland.github.io/pls/reference/mvr.html))

- PLSR - Partial Least Squares Regression
  ([`plsr`](https://khliland.github.io/pls/reference/mvr.html))

- CCA - Canonical Correlation Analysis
  ([`cca`](https://khliland.github.io/multiblock/reference/cca.md))

- IFA - Interbattery Factor Analysis
  ([`ifa`](https://khliland.github.io/multiblock/reference/ifa.md))

- GSVD - Generalized SVD
  ([`gsvd`](https://khliland.github.io/multiblock/reference/gsvd.md))

## See also

Overviews of available methods,
[`multiblock`](https://khliland.github.io/multiblock/reference/multiblock.md),
and methods organised by main structure: `basic`,
[`unsupervised`](https://khliland.github.io/multiblock/reference/unsupervised.md),
[`asca`](https://khliland.github.io/HDANOVA/reference/asca.html),
[`supervised`](https://khliland.github.io/multiblock/reference/supervised.md)
and
[`complex`](https://khliland.github.io/multiblock/reference/complex.md).

## Examples

``` r
data(potato)
X <- potato$Chemical
y <- potato$Sensory[,1,drop=FALSE]

pca.pot  <- pca(X, ncomp = 2)
pcr.pot  <- pcr(y ~ X, ncomp = 2)
pls.pot  <- plsr(y ~ X, ncomp = 2)
cca.pot  <- cca(potato[1:2])
ifa.pot  <- ifa(potato[1:2])
gsvd.pot <- gsvd(lapply(potato[3:4], t))
```
