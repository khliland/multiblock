# Result Functions for Multiblock Objects

Standard result computation and extraction functions for `multiblock`
objects.

## Usage

``` r
# S3 method for class 'multiblock'
scores(object, block = 0, ...)

# S3 method for class 'multiblock'
loadings(object, block = 0, ...)

# S3 method for class 'multiblock'
print(x, ...)

# S3 method for class 'multiblock'
summary(object, ...)
```

## Arguments

- object:

  `multiblock` object.

- block:

  `integer/character` for block selection.

- ...:

  Not implemented.

- x:

  `multiblock` object.

## Value

Scores or loadings are returned by `scores.multiblock` and
`loadings.multiblock`, while print and summary methods invisibly returns
the object.

## Details

Usage of the functions are shown using generics in the examples below.
Object printing and summary are available through: `print.multiblock`
and `summary.multiblock`. Scores and loadings have their own extensions
of [`scores()`](https://khliland.github.io/pls/reference/scores.html)
and [`loadings()`](https://khliland.github.io/pls/reference/scores.html)
throught `scores.multiblock` and `loadings.multiblock`.

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
Common functions for plotting are found in
[`multiblock_plots`](https://khliland.github.io/multiblock/reference/multiblock_plots.md),
respectively.

## Examples

``` r
data(wine)
sc <- sca(wine[c('Smell at rest', 'View', 'Smell after shaking')], ncomp = 4)
print(sc)
#> Simultaneous Component Analysis 
#> 
#> Call:
#> sca(X = wine[c("Smell at rest", "View", "Smell after shaking")],     ncomp = 4)
summary(sc)
#> Simultaneous Component Analysis 
#> =============================== 
#> 
#> $scores: Common scores (21x4)
#> $blockLoadings: Block loadings:
#> - Smell at rest (5x4), View (3x4), Smell after shaking (10x4)
head(loadings(sc, block = 1))
#>                                   Comp 1     Comp 2      Comp 3      Comp 4
#> Odor Intensity before shaking 0.19548267 -0.3267557  0.50138779 -0.09348863
#> Aroma quality before shaking  0.16311889  0.1296282  0.28552072  0.14894779
#> Fruity before shaking         0.12978290  0.1773340  0.34890797  0.26269486
#> Flower before shaking         0.05514227  0.1048857 -0.02986196 -0.47929451
#> Spice before shaking          0.04176668 -0.4175496  0.07332530  0.04104202
head(scores(sc))
#>       Comp 1      Comp 2      Comp 3       Comp 4
#> 1  0.4048505  0.14408310 -0.26713729 -0.221396684
#> 2 -1.2711617  0.25270420 -0.02202010 -0.533496550
#> 3 -0.7036213  0.06526523 -0.11229055  0.078719372
#> 4 -2.0075068 -0.45716224 -0.08211518 -0.007645878
#> 5  1.0467680  0.25398203  0.52254058 -0.029381096
#> 6  0.6670676  0.15451746 -0.58710038  0.050253448
```
