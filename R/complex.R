#' @name complex
#' @title Methods With Complex Linkage
#' @description This documentation covers a few complex methods. In particular:
#' * L-PLS - Partial Least Squares in L configuration (\code{\link{lpls}})
#' * SO-PLS-PM - Sequential and Orthogonalised PLS Path Modeling (\code{\link{sopls_pm}})
#' @examples 
#' # L-PLS
#' sim <- lplsData(I = 30, N = 20, J = 5, K = 6, ncomp = 2)
#' X1  <- sim$X1; X2 <- sim$X2; X3 <- sim$X3
#' lp  <- lpls(X1,X2,X3) # exo-L-PLS
#' 
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
NULL
