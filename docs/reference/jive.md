# Joint and Individual Variation Explained - JIVE

This is a wrapper for the
[`r.jive::jive`](https://rdrr.io/pkg/r.jive/man/jive.html) function for
computing JIVE.

## Usage

``` r
jive(X, ...)
```

## Arguments

- X:

  `list` of input blocks.

- ...:

  additional arguments for
  [`r.jive::jive`](https://rdrr.io/pkg/r.jive/man/jive.html).

## Value

`multiblock` object including relevant scores and loadings. Relevant
plotting functions:
[`multiblock_plots`](https://khliland.github.io/multiblock/reference/multiblock_plots.md)
and result functions:
[`multiblock_results`](https://khliland.github.io/multiblock/reference/multiblock_results.md).

## Details

Jive performs a decomposition of the variation in two or more blocks
into low-dimensional representations of individual and joint variation
plus residual variation.

## References

Lock, E., Hoadley, K., Marron, J., and Nobel, A. (2013) Joint and
individual variation explained (JIVE) for integrated analysis of
multiple data types. Ann Appl Stat, 7 (1), 523–542.

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

## Examples

``` r
 # Too time consuming for testing
  data(candies)
  candyList <- lapply(1:nlevels(candies$candy),function(x)candies$assessment[candies$candy==x,])
  can.jive  <- jive(candyList)
#> Estimating  joint and individual ranks via permutation...
#> Running JIVE algorithm for ranks:
#> joint rank: 1 , individual ranks: 1 1 1 1 1 
#> JIVE algorithm converged after  16  iterations.
#> Re-estimating  joint and individual ranks via permutation...
#> Running JIVE algorithm for ranks:
#> joint rank: 1 , individual ranks: 1 0 1 0 1 
#> JIVE algorithm converged after  13  iterations.
#> Re-estimating  joint and individual ranks via permutation...
#> Final joint rank: 1 , final individual ranks: 1 0 1 0 1 
  summary(can.jive)
#> $Method
#> [1] "perm"
#> 
#> $Ranks
#>      Source     Rank
#> [1,] "Joint"    "1" 
#> [2,] "Source_1" "1" 
#> [3,] "Source_2" "0" 
#> [4,] "Source_3" "1" 
#> [5,] "Source_4" "0" 
#> [6,] "Source_5" "1" 
#> 
#> $Variance
#>            Source_1 Source_2 Source_3 Source_4 Source_5
#> Joint         0.690    0.852    0.793    0.878    0.624
#> Individual    0.079    0.000    0.056    0.000    0.223
#> Residual      0.230    0.148    0.151    0.122    0.153
#> 

```
