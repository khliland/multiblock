# Supervised Multiblock Methods

Collection of supervised multiblock methods:

- MB-PLS - Multiblock Partial Least Squares
  ([`mbpls`](https://khliland.github.io/multiblock/reference/mbpls.md))

- sMB-PLS - Sparse Multiblock Partial Least Squares
  ([`smbpls`](https://khliland.github.io/multiblock/reference/smbpls.md))

- SO-PLS - Sequential and Orthogonalized PLS
  ([`sopls`](https://khliland.github.io/multiblock/reference/sopls.md))

- PO-PLS - Parallel and Orthogonalized PLS
  ([`popls`](https://khliland.github.io/multiblock/reference/popls.md))

- ROSA - Response Oriented Sequential Alternation
  ([`rosa`](https://khliland.github.io/multiblock/reference/rosa.md))

- mbRDA - Multiblock Redundancy Analysis
  ([`mbrda`](https://khliland.github.io/multiblock/reference/mbrda.md))

## See also

Overviews of available methods,
[`multiblock`](https://khliland.github.io/multiblock/reference/multiblock.md),
and methods organised by main structure:
[`basic`](https://khliland.github.io/multiblock/reference/basic.md),
[`unsupervised`](https://khliland.github.io/multiblock/reference/unsupervised.md),
[`asca`](https://khliland.github.io/HDANOVA/reference/asca.html),
`supervised` and
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
mb <- mbpls(Sensory ~ Chemical + Compression, data=potato, ncomp = 5)
print(mb)
#> Multiblock PLS 
#> 
#> Call:
#> mbpls(formula = Sensory ~ Chemical + Compression, data = potato,     ncomp = 5)

# Convert data.frame with AsIs objects to list of matrices
potatoList <- lapply(potato, unclass)
mbr <- mbrda(Sensory ~ Chemical + Compression, data=potatoList, ncomp = 10)
print(mbr)
#> Multiblock RDA 
#> 
#> Call:
#> mbrda(formula = Sensory ~ Chemical + Compression, data = potatoList,     ncomp = 10)
scoreplot(mbr, labels="names")

```
