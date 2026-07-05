# Total, direct, indirect and additional effects in SO-PLS-PM.

SO-PLS-PM is the use of SO-PLS for path-modelling. This particular
function is used to compute effects (explained variances) in sub-paths
of the directed acyclic graph.

## Usage

``` r
sopls_pm(
  X,
  Y,
  ncomp,
  max_comps = min(sum(ncomp), 20),
  sel.comp = "opt",
  computeAdditional = FALSE,
  sequential = FALSE,
  B = NULL,
  k = 10,
  type = "consecutive",
  simultaneous = TRUE
)

# S3 method for class 'SO_TDI'
print(x, showComp = TRUE, heading = "SO-PLS path effects", digits = 2, ...)

sopls_pm_multiple(
  X,
  ncomp,
  max_comps = min(sum(ncomp), 20),
  sel.comp = "opt",
  computeAdditional = FALSE,
  sequential = FALSE,
  B = NULL,
  k = 10,
  type = "consecutive"
)

# S3 method for class 'SO_TDI_multiple'
print(x, heading = "SO-PLS path effects", digits = 2, ...)
```

## Arguments

- X:

  A `list` of input blocks (of type `matrix`).

- Y:

  A `matrix` of response(s).

- ncomp:

  An `integer` vector giving the number of components per block or a
  single integer for common number of components.

- max_comps:

  Maximum total number of components.

- sel.comp:

  A `character` or `integer` vector indicating the type ("opt" - minimum
  error / "chi" - chi-squared reduced) or exact number of components in
  selections.

- computeAdditional:

  A `logical` indicating if additional components should be computed.

- sequential:

  A `logical` indicating if sequential component optimization should be
  applied.

- B:

  An `integer` giving the number of bootstrap replicates for variation
  estimation.

- k:

  An `integer` indicating number of cross-validation segments (default =
  10).

- type:

  A `character` for selecting type of cross-validation segments (default
  = "consecutive").

- simultaneous:

  `logical` indicating if simultaneous orthogonalisation on intermediate
  blocks should be performed (default = TRUE).

- x:

  An object of type `SO_TDI`.

- showComp:

  A `logical` indicating if components should be shown in print (default
  = TRUE).

- heading:

  A `character` giving the heading of the print.

- digits:

  An `integer` for selecting number of digits in print.

- ...:

  Not implemented

## Value

An object of type `SO_TDI` containing total, direct and indirect
effects, plus possibly additional effects and standard deviations
(estimated by bootstrapping).

## Details

`sopls_pm` computes 'total', 'direct', 'indirect' and 'additional'
effects for the 'first' versus the 'last' input block by cross-validated
explained variances. 'total' is the explained variance when doing
regression of 'first' -\> 'last'. 'indirect' is the the same, but
controlled for the intermediate blocks. 'direct' is the left-over part
of the 'total' explained variance when subtracting the 'indirect'.
Finally, 'additional' is the added explained variance of 'last' for each
block following 'first'.

`sopls_pm_multiple` is a wrapper for `sopls_pm` that repeats the
calculation for all pairs of blocks from 'first' to 'last'. Where
`sopls_pm` has a separate response, Y, signifying the 'last' block,
`sopls_pm_multiple` has multiple 'last' blocks, depending on sub-path,
thus collects the response(s) from the list of blocks X.

When sel.comp = "opt", the number of components for all models are
optimized using cross-validation within the ncomp and max_comps
supplied. If sel.comp is "chi", an optimization is also performed, but
parsimonious solutions are sought through a chi-square chriterion. When
setting sel.comp to a numeric vector, exact selection of number of
components is performed.

When setting B to a number, e.g. 200, the procedures above are repeated
B times using bootstrapping to estimate standard deviations of the
cross-validated explained variances.

## References

- Menichelli, E., Almøy, T., Tomic, O., Olsen, N. V., & Næs, T. (2014).
  SO-PLS as an exploratory tool for path modelling. Food quality and
  preference, 36, 122-134.

- Næs, T., Romano, R., Tomic, O., Måge, I., Smilde, A., & Liland, K. H.
  (2020). Sequential and orthogonalized PLS (SO-PLS) regression for path
  analysis: Order of blocks and relations between effects. Journal of
  Chemometrics, e3243.

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
# Single path for the potato data:
data(potato)
pot.pm <- sopls_pm(potato[1:3], potato[['Sensory']], c(5,5,5), computeAdditional=TRUE)
pot.pm
#>  direct indirect     total additional1 additional2 overall
#>   0 (0)    52.44 52.44 (3)    4.09 (3)   14.01 (2)   70.55

# Corresponding SO-PLS model:
# so <- sopls(Sensory ~ ., data=potato[c(1,2,3,9)], ncomp=c(5,5,5), validation="CV", segments=10)
# maageSeq(pot.so, compSeq = c(3,2,4))

# All path in the forward direction for the wine data:
data(wine)
pot.pm.multiple <- sopls_pm_multiple(wine, ncomp = c(4,2,9,8))
pot.pm.multiple
#> $`Smell at rest->View`
#>     direct indirect     total
#>  32.68 (1)        0 32.68 (1)
#> 
#> $`Smell at rest->Smell after shaking`
#>  direct indirect     total
#>   0 (0)    40.03 40.03 (4)
#> 
#> $`Smell at rest->Tasting`
#>  direct indirect     total
#>   0 (0)    11.52 11.52 (2)
#> 
#> $`Smell at rest->Global quality`
#>  direct indirect     total
#>   0 (0)    25.25 25.25 (3)
#> 
#> $`View->Smell after shaking`
#>     direct indirect     total
#>  30.97 (2)        0 30.97 (2)
#> 
#> $`View->Tasting`
#>  direct indirect     total
#>   0 (0)    41.09 41.09 (2)
#> 
#> $`View->Global quality`
#>  direct indirect     total
#>   0 (0)    30.87 30.87 (2)
#> 
#> $`Smell after shaking->Tasting`
#>     direct indirect     total
#>  56.67 (3)        0 56.67 (3)
#> 
#> $`Smell after shaking->Global quality`
#>  direct indirect     total
#>   0 (0)    70.15 70.15 (2)
#> 
#> $`Tasting->Global quality`
#>     direct indirect     total
#>  78.12 (2)        0 78.12 (2)
#> 
```
