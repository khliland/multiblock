# Multiple Factor Analysis - MFA

This is a wrapper for the
[`FactoMineR::MFA`](https://rdrr.io/pkg/FactoMineR/man/MFA.html)
function for computing MFA.

## Usage

``` r
mfa(X, type = rep("c", length(X)), graph = FALSE, ...)
```

## Arguments

- X:

  `list` of input blocks.

- type:

  `character` vector indicating block types, defaults to
  `rep("c", length(X))` for continuous values.

- graph:

  `logical` indicating if decomposition should be plotted.

- ...:

  additional arguments for RGCCA approach.

## Value

`multiblock` object including relevant scores and loadings. Relevant
plotting functions:
[`multiblock_plots`](https://khliland.github.io/multiblock/reference/multiblock_plots.md)
and result functions:
[`multiblock_results`](https://khliland.github.io/multiblock/reference/multiblock_results.md).

## Details

MFA is a methods typically used to compare several equally sized
matrices. It is often used in sensory analyses, where matrices consist
of sensory characteristics and products, and each assessor generates one
matrix each. In its basic form, MFA scales all matrices by their largest
eigenvalue, concatenates them and performs PCA on the result. There are
several possibilities for plots and inspections of the model, handling
of categorical and continuous inputs etc. connected to MFA.

## References

Pagès, J. (2005). Collection and analysis of perceived product
inter-distances using multiple factor analysis: Application to the study
of 10 white wines from the Loire valley. Food Quality and Preference,
16(7), 642–649.

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
pot.mfa    <- mfa(potList)
if(interactive()){
  plot(pot.mfa$MFA)
}
```
