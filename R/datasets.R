#' @name simulated
#' @title Data simulated to have certain characteristics.
#'
#' @description A dataset containing simulated data for 4 connected events where A is the
#' starting point and D is the end point. This can be described as a directed
#' acyclic graph (sketched below, moving left->right). \cr
#' 
#' ![](simulated.png "Path-diagram for simulated data")
#'
#' Subpaths include: ABD, AD, ABCD, ACD
#'
#' @docType data
#' @usage data(simulated)
#' 
#' @format A list of matrices having 200 rows and 10 variables:
#' \describe{
#'   \item{A}{Simulated matrix A}
#'   \item{B}{Simulated matrix B}
#'   ...
#' }
#' 
#' @references Tormod Næs, Rosaria Romano, Oliver Tomic, Ingrid Måge, Age Smilde, Kristian Hovde Liland, 
#' Sequential and orthogonalized PLS (SO-PLS) regression for path analysis: Order of blocks and relations between effects.
#' Journal of Chemometrics, In Press
NULL

#' @name candies
#' @title Sensory assessment of candies.
#'
#' @description A dataset containing 9 sensory attributes for 5 candies assessed
#' by 11 trained assessors.
#'
#' @docType data
#' @usage data(candies)
#' 
#' @format A data.frame having 165 rows and 3 variables:
#' \describe{
#'   \item{assessment}{Matrix of sensory attributes}
#'   \item{assessor}{Factor of assessors}
#'   \item{candy}{Factor of candies}
#' }
#' 
#' @references Luciano G, Næs T. Interpreting sensory data by combining principal 
#' component analysis and analysis of variance. Food Qual Prefer. 2009;20(3):167-175.
NULL

#' @name potato
#' @title Sensory, rheological, chemical and spectroscopic analysis of potatoes.
#'
#' @description A dataset containing 9 blocks of measurements on 26 potatoes.
#' Original dataset can be found at http://models.life.ku.dk/Texture_Potatoes.
#' This version has been pre-processed as follows (corresponding to Liland et al. 2016):
#' * Variables containing NaN have been removed.
#' * Chemical and Compression blocks have been scaled by standard deviations.
#' * NIR blocks have been subjected to SNV (Standard Normal Variate).
#'
#' @docType data
#' @usage data(potato)
#' 
#' @format A data.frame having 26 rows and 9 variables:
#' \describe{
#'   \item{Chemical}{Matrix of chemical measurements}
#'   \item{Compression}{Matrix of rheological compression data}
#'   \item{NIRraw}{Matrix of near-infrared measurements of raw potatoes}
#'   \item{NIRcooked}{Matrix of near-infrared measurements of cooked potatoes}
#'   \item{CPMGraw}{Matrix of NMR (CPMG) measurements of raw potatoes}
#'   \item{CPMGcooked}{Matrix of NMR (CPMG) measurements of cooked potatoes}
#'   \item{FIDraw}{Matrix of NMR (FID) measurements of raw potatoes}
#'   \item{FIDcooked}{Matrix of NMR (FID) measurements of cooked potatoes}
#'   \item{Sensory}{Matrix of sensory assessments}
#' }
#' 
#' @references 
#' * L.G.Thygesen, A.K.Thybo, S.B.Engelsen, Prediction of Sensory Texture Quality of Boiled Potatoes 
#' From Low-field1H NMR of Raw Potatoes. The Role of Chemical Constituents. LWT - Food Science and Technology 34(7), 2001, pp 469-477.
#' * Kristian Hovde Liland, Tormod Næs, Ulf Geir Indahl, ROSA – a fast extension of Partial Least Squares Regression for Multiblock Data Analysis,
#' Journal of Chemometrics 30:11 (2016), pp. 651-662.
NULL

#' @name wine
#' @title Wines of Val de Loire
#'
#' @description This dataset contains sensory assessment of 21 wines. The assessments are grouped
#' according to the tasting process and thus have a natural ordering with a all blocks pointing forward
#' to all remaining blocks in the process.
#' 
#' ![](wine.png "Path-diagram for wine data")
#'
#' @docType data
#' @usage data(wine)
#' 
#' @format A data.frame having 21 rows and 5 variables:
#' \describe{
#'   \item{Smell at rest}{Matrix of sensory assessments}
#'   \item{View}{Matrix of sensory assessments}
#'   \item{Smell after shaking}{Matrix of sensory assessments}
#'   \item{Tasting}{Matrix of sensory assessments}
#'   \item{Global quality}{Matrix of sensory assessments}
#' }
#' 
#' @references Escofier B, Pages L. Analyses Factorielles Simples and Multiples. Paris: Dunod; 1988.
NULL

#' @name mobile
#' @title ECSI Mobile Mobile Phone Provider Dataset
#'
#' @description Mobile data questionnaire often used as an example in path modelling.
#' All the items are scaled from 1 to 10. Score 1 expresses a very negative 
#' point of view on the product while score 10 a very positive opinion. For details,
#' see the original publication.
#' 
#' ![](mobile.png "Path-diagram for mobile data")
#'
#' @docType data
#' @usage data(mobile)
#' 
#' @format A data.frame having 250 rows and 7 variables:
#' \describe{
#'   \item{A}{Image}
#'   \item{B}{Customer expectation}
#'   \item{C}{Perceived quality}
#'   \item{D}{Perceived value}
#'   \item{E}{Customer satisfaction}
#'   \item{F}{Customer complaints}
#'   \item{G}{Customer loyalty}
#' }
#' 
#' @references Tenenhaus M, Esposito Vinzi V, Chatelin YM, Lauro C. PLS path modeling. Comput Stat Data Anal. 2005;48(1):159‐205.
NULL
