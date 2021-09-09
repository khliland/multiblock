#' multiblock
#'
#' @description 
#' A collection of methods for analysis of data sets with more than two blocks of data.
#' 
#' __Unsupervised methods:__
#' * SCA - Simultaneous Component Analysis (\code{\link{sca}})
#' * GCA - Generalized Canonical Analysis (\code{\link{gca}})
#' * GPA - Generalized Procrustes Analysis (\code{\link{gpa}})
#' * MFA - Multiple Factor Analysis (\code{\link{mfa}})
#' * PCA-GCA (\code{\link{pcagca}})
#' * DISCO - Distinctive and Common Components with SCA (\code{\link{disco}})
#' * HPCA - Hierarchical Principal component analysis (\code{\link{hpca}})
#' * MCOA - Multiple Co-Inertia Analysis (\code{\link{mcoa}})
#' * JIVE - Joint and Individual Variation Explained (\code{\link{jive}})
#' * STATIS - Structuration des Tableaux à Trois Indices de la Statistique (\code{\link{statis}})
#' * HOGSVD - Higher Order Generalized SVD (\code{\link{hogsvd}})
#' 
#' __Design based methods:__
#' * ASCA - Anova Simultaneous Component Analysis (\code{\link{asca}})
#' 
#' __Supervised methods:__
#' * MB-PLS - Multiblock Partial Least Squares (\code{\link{mbpls}})
#' * sMB-PLS - Sparse Multiblock Partial Least Squares (\code{\link{smbpls}})
#' * SO-PLS - Sequential and Orthogonalized PLS (\code{\link{sopls}})
#' * PO-PLS - Parallel and Orthogonalized PLS (\code{\link{popls}})
#' * ROSA - Response Oriented Sequential Alternation (\code{\link{rosa}})
#' * mbRDA - Multiblock Redundancy Analysis (\code{\link{mbrda}})
#' 
#' __Complex methods:__
#' * L-PLS - Partial Least Squares in L configuration (\code{\link{lpls}})
#' * SO-PLS-PM - Sequential and Orthogonalised PLS Path Modelling (\code{\link{sopls_pm}})
#'
#' __Single- and two-block methods:__
#' * PCA - Principal Component Analysis (\code{\link{pca}})  
#' * PCR - Principal Component Regression (\code{\link{pcr}})  
#' * PLSR - Partial Least Squares Regression (\code{\link{plsr}})  
#' * CCA - Canonical Correlation Analysis (\code{\link{cca}})  
#' * IFA - Interbattery Factor Analysis (\code{\link{ifa}})
#' * GSVD - Generalized SVD (\code{\link{gsvd}})
#' 
#' __Datasets:__
#' * Sensory assessment of candies (\code{\link{candies}})
#' * Sensory, rheological, chemical and spectroscopic analysis of potatoes (\code{\link{potato}})
#' * Data simulated to have certain characteristics (\code{\link{simulated}})
#' * Wines of Val de Loire (\code{\link{wine}})
#' 
#' __Utility functions:__
#' * Block-wise indexable data.frame (\code{\link{block.data.frame}})
#' * Dummy-code a vector (\code{\link{dummycode}})
#' 
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' @importFrom graphics abline arrows legend
#' @importFrom stats coefficients qf
#' @docType package
#' @name multiblock
NULL