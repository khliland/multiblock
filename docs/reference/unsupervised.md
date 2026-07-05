# Unsupervised Multiblock Methods

Collection of unsupervised multiblock methods:

- SCA - Simultaneous Component Analysis
  ([`sca`](https://khliland.github.io/multiblock/reference/sca.md))

- GCA - Generalized Canonical Analysis
  ([`gca`](https://khliland.github.io/multiblock/reference/gca.md))

- GPA - Generalized Procrustes Analysis
  ([`gpa`](https://khliland.github.io/multiblock/reference/gpa.md))

- MFA - Multiple Factor Analysis
  ([`mfa`](https://khliland.github.io/multiblock/reference/mfa.md))

- PCA-GCA
  ([`pcagca`](https://khliland.github.io/multiblock/reference/pcagca.md))

- DISCO - Distinctive and Common Components with SCA
  ([`disco`](https://khliland.github.io/multiblock/reference/disco.md))

- HPCA - Hierarchical Principal component analysis
  ([`hpca`](https://khliland.github.io/multiblock/reference/hpca.md))

- MCOA - Multiple Co-Inertia Analysis
  ([`mcoa`](https://khliland.github.io/multiblock/reference/mcoa.md))

- JIVE - Joint and Individual Variation Explained
  ([`jive`](https://khliland.github.io/multiblock/reference/jive.md))

- STATIS - Structuration des Tableaux à Trois Indices de la Statistique
  ([`statis`](https://khliland.github.io/multiblock/reference/statis.md))

- HOGSVD - Higher Order Generalized SVD
  ([`hogsvd`](https://khliland.github.io/multiblock/reference/hogsvd.md))

## Details

Original documentation of STATIS:
[statis](https://adeverse.github.io/ade4/reference/statis.html). JIVE,
STATIS and HOGSVD assume variable linked matrices/data.frames, while SCA
handles both links.

## See also

Overviews of available methods,
[`multiblock`](https://khliland.github.io/multiblock/reference/multiblock.md),
and methods organised by main structure:
[`basic`](https://khliland.github.io/multiblock/reference/basic.md),
`unsupervised`,
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
# Object linked data
data(potato)
potList <- as.list(potato[c(1,2,9)])
pot.sca    <- sca(potList)

# Variable linked data
data(candies)
candyList <- lapply(1:nlevels(candies$candy),function(x)candies$assessment[candies$candy==x,])
can.statis <- statis(candyList)
plot(can.statis$statis)
```
