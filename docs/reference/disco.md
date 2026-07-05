# Distinctive and Common Components with SCA - DISCO

This is a wrapper for the `DISCOsca` function by Zhengguo Gu for
computing DISCO.

## Usage

``` r
disco(X, ncomp = 2, ...)
```

## Arguments

- X:

  `list` of input blocks.

- ncomp:

  `integer` number of components to extract.

- ...:

  additional arguments (not used).

## Value

`multiblock` object including relevant scores and loadings. Relevant
plotting functions:
[`multiblock_plots`](https://khliland.github.io/multiblock/reference/multiblock_plots.md)
and result functions:
[`multiblock_results`](https://khliland.github.io/multiblock/reference/multiblock_results.md).

## Details

DISCO is a restriction of SCA where Alternating Least Squares is used
for estimation of loadings and scores. The SCA solution is rotated
towards loadings (in sample linked mode) which are filled with zeros in
a pattern resembling distinct, local and common components. When used in
sample linked mode and only selecting distinct components, it shares a
resemblance to SO-PLS, only in an unsupervised setting. Explained
variances are computed as proportion of block variation explained by
scores\*loadings'.

## References

Schouteden, M., Van Deun, K., Wilderjans, T. F., & Van Mechelen, I.
(2014). Performing DISCO-SCA to search for distinctive and common
information in linked data. Behavior research methods, 46(2), 576-587.

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
data(potato)
potList <- as.list(potato[c(1,2,9)])
pot.disco  <- disco(potList)
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    1    0
#> [3,]    0    1
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
#> [3,]    1    0
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
#> [3,]    0    1
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    0    1
#> [2,]    1    0
#> [3,]    1    0
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    0    1
#> [2,]    1    0
#> [3,]    0    1
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    0    1
#> [2,]    0    1
#> [3,]    1    0
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    1    0
#> [3,]    1    0
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    1    1
#> [3,]    1    0
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    1    0
#> [3,]    1    1
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    1    0
#> [3,]    0    1
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    1    1
#> [3,]    0    1
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    0    1
#> [3,]    1    0
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
#> [3,]    1    1
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    0    1
#> [3,]    0    1
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    0    1
#> [2,]    1    1
#> [3,]    1    0
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    0    1
#> [2,]    1    0
#> [3,]    1    1
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    0    1
#> [2,]    1    1
#> [3,]    0    1
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    0    1
#> [2,]    0    1
#> [3,]    1    1
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    1    1
#> [3,]    1    0
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    1    0
#> [3,]    1    1
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    1    1
#> [3,]    1    1
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    1    1
#> [3,]    0    1
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    0    1
#> [3,]    1    1
#> Now checking the following component structure:
#>      [,1] [,2]
#> [1,]    0    1
#> [2,]    1    1
#> [3,]    1    1
plot(scores(pot.disco), labels="names")

```
