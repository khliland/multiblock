# Sensory, rheological, chemical and spectroscopic analysis of potatoes.

A dataset containing 9 blocks of measurements on 26 potatoes. Original
dataset can be found at http://models.life.ku.dk/Texture_Potatoes. This
version has been pre-processed as follows (corresponding to Liland et
al. 2016):

- Variables containing NaN have been removed.

- Chemical and Compression blocks have been scaled by standard
  deviations.

- NIR blocks have been subjected to SNV (Standard Normal Variate).

## Usage

``` r
data(potato)
```

## Format

A data.frame having 26 rows and 9 variables:

- Chemical:

  Matrix of chemical measurements

- Compression:

  Matrix of rheological compression data

- NIRraw:

  Matrix of near-infrared measurements of raw potatoes

- NIRcooked:

  Matrix of near-infrared measurements of cooked potatoes

- CPMGraw:

  Matrix of NMR (CPMG) measurements of raw potatoes

- CPMGcooked:

  Matrix of NMR (CPMG) measurements of cooked potatoes

- FIDraw:

  Matrix of NMR (FID) measurements of raw potatoes

- FIDcooked:

  Matrix of NMR (FID) measurements of cooked potatoes

- Sensory:

  Matrix of sensory assessments

## References

- L.G.Thygesen, A.K.Thybo, S.B.Engelsen, Prediction of Sensory Texture
  Quality of Boiled Potatoes From Low-field1H NMR of Raw Potatoes. The
  Role of Chemical Constituents. LWT - Food Science and Technology
  34(7), 2001, pp 469-477.

- Kristian Hovde Liland, Tormod Næs, Ulf Geir Indahl, ROSA – a fast
  extension of Partial Least Squares Regression for Multiblock Data
  Analysis, Journal of Chemometrics 30:11 (2016), pp. 651-662.
