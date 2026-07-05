# Inter-battery Factor Analysis - IFA

This is a wrapper for the `RGCCA::rgcca` function for computing IFA.

## Usage

``` r
ifa(X, ncomp = 1, scale = FALSE, verbose = FALSE, ...)
```

## Arguments

- X:

  `list` of input data blocks.

- ncomp:

  `integer` number of principal components to return.

- scale:

  `logical` indicating if variables should be standardised
  (default=FALSE).

- verbose:

  `logical` indicating if intermediate results should be printed.

- ...:

  additional arguments to `RGCCA::rgcca`.

## Value

`multiblock` object with associated with printing, scores, loadings.
Relevant plotting functions:
[`multiblock_plots`](https://khliland.github.io/multiblock/reference/multiblock_plots.md)
and result functions:
[`multiblock_results`](https://khliland.github.io/multiblock/reference/multiblock_results.md).

## Details

IFA rotates two matrices to align one or more factors against each
other, maximising correlations.

## References

Tucker, L. R. (1958). An inter-battery method of factor analysis.
Psychometrika, 23(2), 111-136.

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

ifa.pot  <- ifa(potato[1:2])
#> To run 'ifa', please install the 'RGCCA' package, e.g., using
#> install.packages('RGCCA')
```
