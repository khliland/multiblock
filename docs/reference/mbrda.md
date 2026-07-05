# Multiblock Redundancy Analysis - mbRDA

This is a wrapper for the
[`ade4::mbpcaiv`](https://adeverse.github.io/ade4/reference/mbpcaiv.html)
function for computing mbRDA.

## Usage

``` r
mbrda(formula, data, subset, na.action, X = NULL, Y = NULL, ncomp = 1, ...)
```

## Arguments

- formula:

  Model formula accepting a single response (block) and predictor block
  names separated by + signs.

- data:

  The data set to analyse.

- subset:

  Expression for subsetting the data before modelling.

- na.action:

  How to handle NAs (no action implemented).

- X:

  `list` of input blocks.

- Y:

  `matrix` of responses.

- ncomp:

  `integer` number of PLS components.

- ...:

  additional arguments to ade4::mbpcaiv.

## Value

`multiblock, mvr` object with scores, block-scores and block-loading.
Relevant plotting functions:
[`multiblock_plots`](https://khliland.github.io/multiblock/reference/multiblock_plots.md)
and result functions:
[`multiblock_results`](https://khliland.github.io/multiblock/reference/multiblock_results.md).

## Details

mbRDA is a multiblock formulation of Redundancy (Data) Analysis. RDA is
theoretically between PLS and GCA. Like GCA, RDA does not consider
correlations within X, but like PLS it does consider correlations within
Y. RDA can also be viewed as a PCR of Y constrained to have scores that
are also linear combinations of X. If the `adegraphics` package is
attached, a nice overview can be plotted as `plot(mbr$mbpcaiv)`
following the example below.

## References

Bougeard, S., Qannari, E.M., Lupo, C., andHanafi, M. (2011). From
Multiblock Partial Least Squares to Multiblock Redundancy Analysis. A
Continuum Approach. Informatica, 22(1), 11–26.

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
# Convert data.frame with AsIs objects to list of matrices
data(potato)
potatoList <- lapply(potato, unclass)

mbr <- mbrda(Sensory ~ Chemical + Compression, data = potatoList, ncomp = 10)
mbr.XY <- mbrda(X = potatoList[c('Chemical','Compression')], Y = potatoList[['Sensory']], 
                ncomp = 10)
print(mbr)
#> Multiblock RDA 
#> 
#> Call:
#> mbrda(formula = Sensory ~ Chemical + Compression, data = potatoList,     ncomp = 10)
scoreplot(mbr) # Exploiting mvr object structure from pls package
```
