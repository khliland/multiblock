# Block-wise indexable data.frame

This is a convenience function for making `data.frame`s that are easily
indexed on a block-wise basis.

## Usage

``` r
block.data.frame(X, block_inds = NULL, to.matrix = TRUE)
```

## Arguments

- X:

  Either a single `data.frame` to index or a `list` of
  matrices/data.frames

- block_inds:

  Named `list` of indexes if `X` is a single `data.frame`, otherwise
  `NULL`.

- to.matrix:

  `logical` indicating if input list elements should be converted to
  matrices.

## Value

A `data.frame` which can be indexed block-wise.

## Examples

``` r
# Random data
M <- matrix(rnorm(200), nrow = 10)
# .. with dimnames
dimnames(M) <- list(LETTERS[1:10], as.character(1:20))

# A named list for indexing
inds <- list(B1 = 1:10, B2 = 11:20)

X <- block.data.frame(M, inds)
str(X)
#> 'data.frame':    10 obs. of  2 variables:
#>  $ B1: 'AsIs' num [1:10, 1:10] -0.936 -0.016 -0.827 -1.512 0.935 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:10] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:10] "1" "2" "3" "4" ...
#>  $ B2: 'AsIs' num [1:10, 1:10] 0.604 -0.263 -0.528 0.192 -1.146 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:10] "A" "B" "C" "D" ...
#>   .. ..$ : chr [1:10] "11" "12" "13" "14" ...
```
