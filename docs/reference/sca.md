# Simultaneous Component Analysis - SCA

This is a basic implementation of the SCA-P algorithm (least restricted
SCA) with support for both sample- and variable-linked modes.

## Usage

``` r
sca(X, ncomp = 2, scale = FALSE, samplelinked = "auto", ...)
```

## Arguments

- X:

  `list` of input blocks.

- ncomp:

  `integer` number of components to extract.

- scale:

  `logical` indicating autoscaling of features (default = FALSE).

- samplelinked:

  `character/logical` indicating if blocks are linked by samples (TRUE)
  or variables (FALSE). Using 'auto' (default), this will be determined
  automatically.

- ...:

  additional arguments (not used).

## Value

`multiblock` object including relevant scores and loadings. Relevant
plotting functions:
[`multiblock_plots`](https://khliland.github.io/multiblock/reference/multiblock_plots.md)
and result functions:
[`multiblock_results`](https://khliland.github.io/multiblock/reference/multiblock_results.md).

## Details

SCA, in its original variable-linked version, calculates common loadings
and block-wise scores. There are many possible constraints and
specialisations. This implementations uses PCA as the backbone, thus
resulting in deterministic, ordered components. A parameter controls the
linking mode, but if left untouched an attempt is made at automatically
determining variable or sample linking.

## References

Levin, J. (1966) Simultaneous factor analysis of several gramian
matrices. Psychometrika, 31(3), 413–419.

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
# Object linked data
data(potato)
potList <- as.list(potato[c(1,2,9)])
pot.sca    <- sca(potList)
plot(scores(pot.sca), labels="names")


# Variable linked data
data(candies)
candyList <- lapply(1:nlevels(candies$candy),function(x)candies$assessment[candies$candy==x,])
pot.sca    <- sca(candyList, samplelinked = FALSE)
pot.sca
#> Simultaneous Component Analysis 
#> 
#> Call:
#> sca(X = candyList, samplelinked = FALSE)
```
