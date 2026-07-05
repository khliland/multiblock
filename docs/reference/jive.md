# Joint and Individual Variation Explained - JIVE

This is a wrapper for the `r.jive::jive` function for computing JIVE.

## Usage

``` r
jive(X, ...)
```

## Arguments

- X:

  `list` of input blocks.

- ...:

  additional arguments for `r.jive::jive`.

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
#> To run 'jive', please install the 'r.jive' package, e.g., using
#> install.packages('r.jive')
  summary(can.jive)
#> Nothing computed 
#> ================ 
#> 
#> $scores: Not used (2x2)
#> $loadings: Not used (2x2)
#> $blockScores: Not used:
#> - (2x2), (2x2), (2x2), (2x2), (2x2)

```
