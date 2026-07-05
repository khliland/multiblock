# A. Data handling

``` r

# Start the multiblock R package
library(multiblock)
#> Registered S3 method overwritten by 'lme4':
#>   method           from
#>   na.action.merMod car
#> Registered S3 method overwritten by 'plsVarSel':
#>   method       from
#>   print.mvrVal pls
#> Registered S3 methods overwritten by 'multiblock':
#>   method             from
#>   print.multiblock   ade4
#>   summary.multiblock ade4
#> 
#> Attaching package: 'multiblock'
#> The following object is masked from 'package:stats':
#> 
#>     loadings
```

## Read from file

Data are stored in many different file formats. The following three
examples cover two types of CSV-files and generic flat files.

``` r

# Find directory extdata from the multiblock package
mbdir <- system.file('extdata/', package = "multiblock")

# Comma separated values, row names in first column
meta_data <- read.csv(paste0(mbdir, "/meta_data.csv"), row.names = 1)
# If working directory matches file location:
# meta_data <- read.csv('meta_data.csv', row.names = 1)
meta_data
#>         temperature colour
#> John           38.0   blue
#> Julia          37.0  green
#> James          37.5   blue
#> Jacob          37.6    red
#> Jane           37.2    red
#> Johanna        37.9  green

# Semi-colon separated values (locales where the decimal point is comma),
# no row names
proteins <- read.csv2(paste0(mbdir, "/proteins.csv"))
proteins
#>         prot1       prot2      prot3
#> 1  0.46532048  0.30183300 -1.4654414
#> 2 -1.79802081 -0.22812232 -0.4639203
#> 3 -1.92962434 -0.40513080  0.1767796
#> 4  0.87437138  0.79843798  0.1234731
#> 5 -0.62445278 -0.07975479 -1.1126332
#> 6 -0.07493721  1.09576027  1.2656596

# Blank space separated data without labels
genes <- read.table(paste0(mbdir, "/genes.dat"))
genes
#>            V1         V2         V3
#> 1  0.39033106 -0.5720390  1.9147573
#> 2  0.55352785  0.0948703 -0.2239755
#> 3  0.09872346 -0.1029385  0.9047138
#> 4 -0.59213740 -0.6027739  0.6177083
#> 5 -0.02350148  0.3572809 -0.5168416
#> 6  0.76644845  1.2863428  1.8239298
```

## Data pre-processing

Before analysis, various types of pre-processing may be needed. Centring
and standardising/scaling may be considered the most basic. In R, these
operations are performed column-wise by default, leading to autoscaling.
If these operations are performed on the rows, we perform the standard
normal variate (SNV) instead.

``` r

# Column-centring
genes_centred <- scale(genes, scale=FALSE)
colMeans(genes_centred) # Check mean values
#>           V1           V2           V3 
#> 1.850372e-17 0.000000e+00 7.401487e-17

# Autoscaling
genes_scaled <- scale(genes)
apply(genes_scaled, 2, sd) # Check standard deviations
#> V1 V2 V3 
#>  1  1  1

# SNV (transpose, autoscale, re-transpose)
genes_snv <- t(scale(t(genes)))
apply(genes_snv, 1, sd) # Check standard deviations
#> [1] 1 1 1 1 1 1
```

### Re-coding categorical data

Most analysis methods require continuous input data. The **meta_data**
**data.frame** contains a character vector (a factor in older R
versions) of categories. This package has a function **dummycode** for
converting categorical data to various dummy formats.

``` r

# Default is sum coding
dummycode(meta_data$colour)
#>   x1 x2
#> 1  1  0
#> 2  0  1
#> 3  1  0
#> 4 -1 -1
#> 5 -1 -1
#> 6  0  1

# Treatment coding
dummycode(meta_data$colour, "contr.treatment")
#>   xgreen xred
#> 1      0    0
#> 2      1    0
#> 3      0    0
#> 4      0    1
#> 5      0    1
#> 6      1    0

# Full dummy-coding (rank deficient)
dummycode(meta_data$colour, drop = FALSE)
#>   xblue xgreen xred
#> 1     1      0    0
#> 2     0      1    0
#> 3     1      0    0
#> 4     0      0    1
#> 5     0      0    1
#> 6     0      1    0

# Replace categorical with dummy-coded, use I() to index by common name
meta_data2 <- meta_data
meta_data2$colour <- I(dummycode(meta_data$colour, drop = FALSE))
meta_data2
#>         temperature colour.xblue colour.xgreen colour.xred
#> John           38.0            1             0           0
#> Julia          37.0            0             1           0
#> James          37.5            1             0           0
#> Jacob          37.6            0             0           1
#> Jane           37.2            0             0           1
#> Johanna        37.9            0             1           0
meta_data2$colour
#>   xblue xgreen xred
#> 1     1      0    0
#> 2     0      1    0
#> 3     1      0    0
#> 4     0      0    1
#> 5     0      0    1
#> 6     0      1    0
```

## Data structures for multiblock analysis

### Create list of blocks

A simple list of blocks can be created using the **list()** function.
Naming of the blocks can be done directly or after creation.

``` r

# Direct approach
blocks1 <- list(meta = meta_data2, proteins = proteins, genes = genes)

# Two-step approach
blocks2 <- list(meta_data2, proteins, genes)
names(blocks2) <- c('meta', 'proteins', 'genes')

# Same result
identical(blocks1, blocks2)
#> [1] TRUE

# Access by name or number
blocks1[['meta']]
#>         temperature colour.xblue colour.xgreen colour.xred
#> John           38.0            1             0           0
#> Julia          37.0            0             1           0
#> James          37.5            1             0           0
#> Jacob          37.6            0             0           1
#> Jane           37.2            0             0           1
#> Johanna        37.9            0             1           0
blocks2[[1]]
#>         temperature colour.xblue colour.xgreen colour.xred
#> John           38.0            1             0           0
#> Julia          37.0            0             1           0
#> James          37.5            1             0           0
#> Jacob          37.6            0             0           1
#> Jane           37.2            0             0           1
#> Johanna        37.9            0             1           0
```

### Create data.frame of blocks

A **data.frame** is a convenient storage format for data in R and can
handle many types of variables, e.g. numeric, logical, character, factor
or matrices. The latter is useful for analyses of data with shared
sample mode.

``` r

# Construct block data.frame from list
df1 <- block.data.frame(blocks1)

# Construct block data.frame from data.frame:
# First merge blocks into data.frame
my_data <- cbind(meta_data2, proteins, genes)
# Then construct block data.frame using named 
# list of indexes
df2 <- block.data.frame(my_data, block_inds = 
        list(meta = 1:2, proteins = 3:5, genes = 6:8))

# Same result
identical(df1,df2)
#> [1] TRUE

# Access by name or number
df1[[2]]
#>               prot1       prot2      prot3
#> John     0.46532048  0.30183300 -1.4654414
#> Julia   -1.79802081 -0.22812232 -0.4639203
#> James   -1.92962434 -0.40513080  0.1767796
#> Jacob    0.87437138  0.79843798  0.1234731
#> Jane    -0.62445278 -0.07975479 -1.1126332
#> Johanna -0.07493721  1.09576027  1.2656596
df2[['proteins']]
#>               prot1       prot2      prot3
#> John     0.46532048  0.30183300 -1.4654414
#> Julia   -1.79802081 -0.22812232 -0.4639203
#> James   -1.92962434 -0.40513080  0.1767796
#> Jacob    0.87437138  0.79843798  0.1234731
#> Jane    -0.62445278 -0.07975479 -1.1126332
#> Johanna -0.07493721  1.09576027  1.2656596
df1[c(1,3)]
#>         meta.temperature meta.colour.xblue meta.colour.xgreen meta.colour.xred
#> John                38.0               1.0                0.0              0.0
#> Julia               37.0               0.0                1.0              0.0
#> James               37.5               1.0                0.0              0.0
#> Jacob               37.6               0.0                0.0              1.0
#> Jane                37.2               0.0                0.0              1.0
#> Johanna             37.9               0.0                1.0              0.0
#>            genes.V1    genes.V2    genes.V3
#> John     0.39033106 -0.57203901  1.91475732
#> Julia    0.55352785  0.09487030 -0.22397553
#> James    0.09872346 -0.10293847  0.90471382
#> Jacob   -0.59213740 -0.60277392  0.61770832
#> Jane    -0.02350148  0.35728086 -0.51684156
#> Johanna  0.76644845  1.28634283  1.82392983
df1[-2]
#>         meta.temperature meta.colour.xblue meta.colour.xgreen meta.colour.xred
#> John                38.0               1.0                0.0              0.0
#> Julia               37.0               0.0                1.0              0.0
#> James               37.5               1.0                0.0              0.0
#> Jacob               37.6               0.0                0.0              1.0
#> Jane                37.2               0.0                0.0              1.0
#> Johanna             37.9               0.0                1.0              0.0
#>            genes.V1    genes.V2    genes.V3
#> John     0.39033106 -0.57203901  1.91475732
#> Julia    0.55352785  0.09487030 -0.22397553
#> James    0.09872346 -0.10293847  0.90471382
#> Jacob   -0.59213740 -0.60277392  0.61770832
#> Jane    -0.02350148  0.35728086 -0.51684156
#> Johanna  0.76644845  1.28634283  1.82392983
df2[c('proteins','genes')]
#>         proteins.prot1 proteins.prot2 proteins.prot3    genes.V1    genes.V2
#> John        0.46532048     0.30183300    -1.46544136  0.39033106 -0.57203901
#> Julia      -1.79802081    -0.22812232    -0.46392027  0.55352785  0.09487030
#> James      -1.92962434    -0.40513080     0.17677963  0.09872346 -0.10293847
#> Jacob       0.87437138     0.79843798     0.12347308 -0.59213740 -0.60277392
#> Jane       -0.62445278    -0.07975479    -1.11263323 -0.02350148  0.35728086
#> Johanna    -0.07493721     1.09576027     1.26565956  0.76644845  1.28634283
#>            genes.V3
#> John     1.91475732
#> Julia   -0.22397553
#> James    0.90471382
#> Jacob    0.61770832
#> Jane    -0.51684156
#> Johanna  1.82392983

# Use with formula interface (see other vignettes)
# sopls(meta ~ proteins + genes, data = df1)

# Use with single list interface (see other vignettes)
# mfa(df1[c(1,3)], ncomp = 3)
```
