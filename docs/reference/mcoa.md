# Multiple Co-Inertia Analysis - MCOA

This is a wrapper for the `RGCCA::rgcca` function for computing MCOA.

## Usage

``` r
mcoa(X, ncomp = 2, scale = FALSE, verbose = FALSE, ...)
```

## Arguments

- X:

  `list` of input blocks.

- ncomp:

  `integer` number of components to extract.

- scale:

  `logical` indicating if variables should be scaled.

- verbose:

  `logical` indicating if diagnostic information should be printed.

- ...:

  additional arguments for RGCCA.

## Value

`multiblock` object including relevant scores and loadings. Relevant
plotting functions:
[`multiblock_plots`](https://khliland.github.io/multiblock/reference/multiblock_plots.md)
and result functions:
[`multiblock_results`](https://khliland.github.io/multiblock/reference/multiblock_results.md).

## Details

MCOA resembles GCA and MFA in that it creates a set of reference scores,
for which each block's individual scores should correlate maximally too,
but also the variance within each block should be taken into account. A
single component solution is equivalent to a PCA on concatenated blocks
scaled by the so called inverse inertia.

## References

- Le Roux; B. and H. Rouanet (2004). Geometric Data Analysis, From
  Correspondence Analysis to Structured Data Analysis. Dordrecht.
  Kluwer: p.180.

- Greenacre, Michael and Blasius, Jörg (editors) (2006). Multiple
  Correspondence Analysis and Related Methods. London: Chapman &
  Hall/CRC.

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
pot.mcoa   <- mcoa(potList)
#> To run 'mcoa', please install the 'RGCCA' package, e.g., using
#> install.packages('RGCCA')
plot(scores(pot.mcoa), labels="names")

```
