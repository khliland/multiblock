#' multiblock
#'
#' 
#' @description 
#' A collection of methods for analysis of data sets with more than two blocks of data.
#' 
#' __Unsupervised methods:__
#' * SCA - Simultaneous Component Analysis (\code{sca})
#' * GCA - Generalized Canonical Analysis (\code{gca})
#' * GPA - Generalized Procrustes Analysis (\code{gpa})
#' * MFA - Multiple Factor Analysis (\code{mfa})
#' * PCA-GCA (\code{pcagca})
#' * DISCO - Distinctive and Common Components with SCA (\code{disco})
#' * HPCA - Hierarchical Principal component analysis (\code{hpca})
#' * MCOA - Multiple Co-Inertia Analysis (\code{mcoa})
#' * JIVE - Joint and Individual Variation Explained (\code{jive})
#' * STATIS - Structuration des Tableaux à Trois Indices de la Statistique (\code{statis})
#' * HOGSVD - Higher Order Generalized SVD (\code{hogsvd})
#' 
#' __Design based methods:__
#' * ASCA - Anova Simultaneous Component Analysis (\code{asca})
#' 
#' __Supervised methods:__
#' * MB-PLS - Multiblock Partial Least Squares (\code{mbpls})
#' * SO-PLS - Sequential and Orthogonalized PLS (\code{\link{sopls}})
#' * PO-PLS - Parallel and Orthogonalized PLS (\code{popls})
#' * ROSA - Response Oriented Sequential Alternation (\code{\link{rosa}})
#' * mbRDA - Multiblock Redundancy Analysis (\code{mbrda})
#' 
#' __Complex methods:__
#' * L-PLS - Partial Least Squares in L configuration (_lpls_)
#' * SO-PLS-PM - Sequential and Orthogonalised PLS Path Modeling (_sopls_pm_)
#'
#' __Single- and two-block methods:__
#' * PCA - Principal Component Analysis (\code{pca})  
#' * PCR - Principal Component Regression (\code{pcr})  
#' * PLSR - Partial Least Squares Regression (\code{plsr})  
#' * CCA - Canonical Correlation Analysis (\code{cca})  
#' * IFA - Interbattery Factor Analysis (\code{ifa})
#' * GSVD - Generalized SVD (\code{gsvd})
#' 
#' @seealso Overviews of available methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' @importFrom graphics abline arrows legend
#' @importFrom stats coefficients qf
#' @docType package
#' @name multiblock
NULL