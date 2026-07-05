# Response Oriented Sequential Alternation - ROSA

Formula based interface to the ROSA algorithm following the style of the
`pls` package.

## Usage

``` r
rosa(
  formula,
  ncomp,
  Y.add,
  common.comp = 1,
  data,
  subset,
  na.action,
  scale = FALSE,
  weights = NULL,
  validation = c("none", "CV", "LOO"),
  internal.validation = FALSE,
  fixed.block = NULL,
  design.block = NULL,
  canonical = TRUE,
  ...
)
```

## Arguments

- formula:

  Model formula accepting a single response (block) and predictor block
  names separated by + signs.

- ncomp:

  The maximum number of ROSA components.

- Y.add:

  Optional response(s) available in the data set.

- common.comp:

  Automatically create all combinations of common components up to
  length `common.comp` (default = 1).

- data:

  The data set to analyse.

- subset:

  Expression for subsetting the data before modelling.

- na.action:

  How to handle NAs (no action implemented).

- scale:

  Optionally scale predictor variables by their individual standard
  deviations.

- weights:

  Optional object weights.

- validation:

  Optional cross-validation strategy "CV" or "LOO".

- internal.validation:

  Optional cross-validation for block selection process, "LOO", "CV3",
  "CV5", "CV10" (CV-number of segments), or vector of integers (default
  = FALSE).

- fixed.block:

  integer vector with block numbers for each component (0 = not fixed)
  or list of length \<= ncomp (element length 0 = not fixed).

- design.block:

  integer vector containing block numbers of design blocks

- canonical:

  logical indicating if canonical correlation should be use when
  calculating loading weights (default), enabling B/W maximization,
  common components, etc. Alternatively (FALSE) a PLS2 strategy, e.g.
  for spectra response, is used.

- ...:

  Additional arguments for `cvseg` or `rosa.fit`

## Value

An object of classes `rosa` and `mvr` having several associated printing
([`rosa_results`](https://khliland.github.io/multiblock/reference/rosa_results.md))
and plotting methods
([`rosa_plots`](https://khliland.github.io/multiblock/reference/rosa_plots.md)).

## Details

ROSA is an opportunistic method sequentially selecting components from
whichever block explains the response most effectively. It can be
formulated as a PLS model on concatenated input block with block
selection per component. This implementation adds several options that
are not described in the literature. Most importantly, it opens for
internal validation in the block selection process, making this more
robust. In addition it handles design blocks explicitly, enables
classification and secondary responses (CPLS), and definition of common
components.

## References

Liland, K.H., Næs, T., and Indahl, U.G. (2016). ROSA - a fast extension
of partial least squares regression for multiblock data analysis.
Journal of Chemometrics, 30, 651–662, doi:10.1002/cem.2824.

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
[`rosa_results`](https://khliland.github.io/multiblock/reference/rosa_results.md)
and
[`rosa_plots`](https://khliland.github.io/multiblock/reference/rosa_plots.md),
respectively.

## Examples

``` r
data(potato)
mod <- rosa(Sensory[,1] ~ ., data = potato, ncomp = 10, validation = "CV", segments = 5)
summary(mod)
#> Data:    X dimension: 26 3946 
#>  Y dimension: 26 1
#> Fit method:
#> Number of components considered: 10
#> 
#> VALIDATION: RMSEP
#> Cross-validated using 5 random segments.
#>        (Intercept)  1 comps  2 comps  3 comps  4 comps  5 comps  6 comps
#> CV           1.778    1.421    1.127   0.8676   0.7277   0.6523   0.5383
#> adjCV        1.778    1.338    1.061   0.7969   0.6612   0.5926   0.4969
#>        7 comps  8 comps  9 comps  10 comps
#> CV      0.4826   0.4652   0.4928    0.5961
#> adjCV   0.4482   0.4303   0.4497    0.5483
#> 
#> TRAINING: % variance explained
#>               1 comps  2 comps  3 comps  4 comps  5 comps  6 comps  7 comps
#> X               27.82    40.72    47.95    49.48    53.96    55.28    57.79
#> Sensory[, 1]    76.54    86.77    94.07    96.47    97.03    97.49    97.83
#>               8 comps  9 comps  10 comps
#> X               65.49    67.14     68.16
#> Sensory[, 1]    98.03    98.28     98.51

# For examples of ROSA results and plotting see 
# ?rosa_results and ?rosa_plots.
```
