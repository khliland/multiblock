# Unique combinations of blocks

Compute a list of all possible block combinations where the number of
blocks in each combination is limited by parameters `min_level` and
`max_level`.

## Usage

``` r
unique_combos(n_block, max_level, min_level = 2)
```

## Arguments

- n_block:

  `integer` number of input blocks.

- max_level:

  `integer` maximum number of blocks per combination.

- min_level:

  `integer` minimum number of blocks per combination.

## Value

A list of unique block combinations.

## Details

This function is used for minimal redundancy implementations of
[`rosa`](https://khliland.github.io/multiblock/reference/rosa.md) and
[`sopls`](https://khliland.github.io/multiblock/reference/sopls.md) and
for lookups into computed components.

## Examples

``` r
unique_combos(3, 2)
#> [[1]]
#> [1] 1 2
#> 
#> [[2]]
#> [1] 1 3
#> 
#> [[3]]
#> [1] 2 3
#> 
```
