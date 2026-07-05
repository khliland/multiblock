# Dummy-coding of a single vector

Flexible dummy-coding allowing for all R's built-in types of contrasts
and optional dropping of a factor level to reduce rank defficiency
probability.

## Usage

``` r
dummycode(Y, contrast = "contr.sum", drop = TRUE)
```

## Arguments

- Y:

  `vector` to dummy code.

- contrast:

  Contrast type, default = "contr.sum".

- drop:

  `logical` indicating if one level should be dropped (default = TRUE).

## Value

`matrix` made by dummy-coding the input vector.

## Examples

``` r
vec <- c("a","a","b","b","c","c")
dummycode(vec)
#>   x1 x2
#> 1  1  0
#> 2  1  0
#> 3  0  1
#> 4  0  1
#> 5 -1 -1
#> 6 -1 -1
```
