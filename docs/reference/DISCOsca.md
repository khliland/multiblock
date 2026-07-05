# DISCO-SCA rotation.

A DISCO-SCA procedure for identifying common and distinctive components.
The code is adapted from the orphaned RegularizedSCA package by Zhengguo
Gu.

## Usage

``` r
DISCOsca(DATA, R, Jk)
```

## Arguments

- DATA:

  A matrix, which contains the concatenated data with the same subjects
  from multiple blocks. Note that each row represents a subject.

- R:

  Number of components (R\>=2).

- Jk:

  A vector containing number of variables in the concatenated data
  matrix.

## Value

- Trot_best:

  Estimated component score matrix (i.e., T)

- Prot_best:

  Estimated component loading matrix (i.e., P)

- comdist:

  A matrix representing common distinctive components. (Rows are data
  blocks and columns are components.) 0 in the matrix indicating that
  the corresponding component of that block is estimated to be zeros,
  and 1 indicates that (at least one component loading in) the
  corresponding component of that block is not zero. Thus, if a column
  in the `comdist` matrix contains only 1's, then this column is a
  common component, otherwise distinctive component.

- propExp_component:

  Proportion of variance per component.

## References

Schouteden, M., Van Deun, K., Wilderjans, T. F., & Van Mechelen, I.
(2014). Performing DISCO-SCA to search for distinctive and common
information in linked data. Behavior research methods, 46(2), 576-587.

## Examples

``` r
if (FALSE) { # \dontrun{
DATA1 <- matrix(rnorm(50), nrow=5)
DATA2 <- matrix(rnorm(100), nrow=5) 
DATA <- cbind(DATA1, DATA2)
R <- 5
Jk <- c(10, 20) 
DISCOsca(DATA, R, Jk)
} # }
```
