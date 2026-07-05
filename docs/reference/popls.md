# Parallel and Orthogonalised Partial Least Squares - PO-PLS

This is a basic implementation of PO-PLS with manual and automatic
component selections.

## Usage

``` r
popls(
  X,
  Y,
  commons = 2,
  auto = TRUE,
  auto.par = list(explVarLim = 40, rLim = 0.8),
  manual.par = list(ncomp = rep(0, length(X)), ncommon = list())
)
```

## Arguments

- X:

  `list` of input blocks

- Y:

  `matrix` of response variable(s)

- commons:

  `numeric` giving the highest number of blocks to combine when
  calculating local or common scores.

- auto:

  `logical` indicating if automatic choice of complexities should be
  used.

- auto.par:

  `named list` setting limits for automatic choice of complexities. See
  Details.

- manual.par:

  `named list` for manual choice of blocks. The list consists of `ncomp`
  which indicates the number of components to extract from each block
  and `ncommon` which is the corresponding for choosing the block
  combinations (local/common). For the latter, use
  unique_combos(n_blocks, commons) to see order of local/common blocks.
  Component numbers will be reduced if simpler models give better
  predictions. See example.

## Value

A `multiblock` object with block-wise, local and common loadings and
scores. Relevant plotting functions:
[`multiblock_plots`](https://khliland.github.io/multiblock/reference/multiblock_plots.md)
and result functions:
[`multiblock_results`](https://khliland.github.io/multiblock/reference/multiblock_results.md).

## Details

PO-PLS decomposes a set of input data blocks into common, local and
distinct components through a process involving `pls` and
[`gca`](https://khliland.github.io/multiblock/reference/gca.md). The
`rLim` parameter is a lower bound for the GCA correlation when building
common components, while explVarLim is the minimum explained variance
for common components and unique components.

## References

- I Måge, BH Mevik, T Næs. (2008). Regression models with process
  variables and parallel blocks of raw material measurements. Journal of
  Chemometrics: A Journal of the Chemometrics Society 22 (8), 443-456

- I Måge, E Menichelli, T Næs (2012). Preference mapping by PO-PLS:
  Separating common and unique information in several data blocks. Food
  quality and preference 24 (1), 8-16

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

# Automatic analysis
pot.po.auto <- popls(potato[1:3], potato[['Sensory']][,1])
#> Warning: 'ncomp' reduced due to low singular value for block 1
#> Warning: 'ncomp' reduced due to low singular value for block 1
pot.po.auto$explVar
#> $Chemical
#> named numeric(0)
#> 
#> $Compression
#> C(1,2), Comp 1 
#>       68.14595 
#> 
#> $NIRraw
#> D(3), Comp 1 
#>     42.53836 
#> 

# Manual choice of up to 5 components for each block and 1, 0, and 2 blocks,
# respectively from the (1,2), (1,3) and (2,3) combinations of blocks.
pot.po.man <- popls(potato[1:3], potato[['Sensory']][,1], auto=FALSE, 
                manual.par = list(ncomp=c(5,5,5), ncommon=c(1,0,2)))
#> Warning: 'ncomp' reduced due to low singular value for block 1
#> Warning: 'ncomp' reduced due to low singular value for block 1
pot.po.man$explVar
#> $Chemical
#> C(1,2), Comp 1   D(1), Comp 1   D(1), Comp 2   D(1), Comp 3 
#>    32.22943509    20.74356738    20.87956167     0.04446336 
#> 
#> $Compression
#> C(1,2), Comp 1 C(2,3), Comp 1 C(2,3), Comp 2 
#>      68.145954       5.504006       7.214911 
#> 
#> $NIRraw
#> C(2,3), Comp 1 C(2,3), Comp 2 
#>      32.606459       5.475702 
#> 

# Score plot for local (2,3) components
plot(scores(pot.po.man,3), comps=1:2, labels="names")

```
