# PCA-GCA

PCA-GCA is a methods which aims at estimating subspaces of common, local
and distinct variation from two or more blocks.

## Usage

``` r
pcagca(
  X,
  commons = 2,
  auto = TRUE,
  auto.par = list(explVarLim = 40, rLim = 0.8),
  manual.par = list(ncomp = 0, ncommon = 0),
  tol = 10^-12
)
```

## Arguments

- X:

  `list` of input blocks

- commons:

  `numeric` giving the highest number of blocks to combine when
  calculating local or common scores.

- auto:

  `logical` indicating if automatic choice of complexities should be
  used.

- auto.par:

  `named list` setting limits for automatic choice of complexities.

- manual.par:

  `named list` for manual choice of blocks. The list consists of `ncomp`
  which indicates the number of components to extract from each block
  and `ncommon` which is the corresponding for choosing the block
  combinations (local/common). For the latter, use
  unique_combos(n_blocks, commons) to see order of local/common blocks.
  Component numbers will be reduced if simpler models give better
  predictions. See example.

- tol:

  `numeric` tolerance for component inclusion (singular values).

## Value

`multiblock` object including relevant scores and loadings. Relevant
plotting functions:
[`multiblock_plots`](https://khliland.github.io/multiblock/reference/multiblock_plots.md)
and result functions:
[`multiblock_results`](https://khliland.github.io/multiblock/reference/multiblock_results.md).
Distinct components are marked as 'D(x), Comp c' for block x and
component c while local and common components are marked as "C(x1, x2),
Comp c", where x1 and x2 (and more) are block numbers.

## Details

The name PCA-GCA comes from the process of first applying PCA to each
block, then using GCA to estimate local and common components, and
finally orthogonalising the block-wise scores on the local/common ones
and re-estimating these to obtain distinct components. The procedure is
highly similar to the supervised method PO-PLS, where the PCA steps are
exchanged with PLS.

## References

Smilde, A., Måge, I., Naes, T., Hankemeier, T.,Lips, M., Kiers, H.,
Acar, E., and Bro, R.(2017). Common and distinct components in data
fusion. Journal of Chemometrics, 31(7), e2900.

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
pot.pcagca <- pcagca(potList)
#> Warning: 'ncomp' reduced due to low singular value for block 2

# Show origin and type of all components
lapply(pot.pcagca$blockScores,colnames)
#> $Chemical
#> NULL
#> 
#> $Compression
#> [1] "D(2), Comp 1"
#> 
#> $Sensory
#> [1] "C(1,3), Comp 1"
#> 

# Basic multiblock plot
plot(scores(pot.pcagca, block=2), comps=1, labels="names")

```
