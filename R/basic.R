#' @name basic
#' @title Single- and Two-Block Methods
#' @description This documentation covers a range of single- and two-block methods. In particular:
#' * PCA - Principal Component Analysis (\code{pca})  
#' * PCR - Principal Component Regression (\code{pcr})  
#' * PLSR - Partial Least Squares Regression (\code{plsr})  
#' * CCA - Canonical Correlation Analysis (\code{cca})  
#' * IFA - Interbattery Factor Analysis (\code{ifa})
#' * GSVD - Generalized SVD (\code{gsvd})
#' 
#' @importFrom pls pcr plsr
#' @importFrom geigen gsvd
#' @rdname basic 
#' 
#' @references 
#' * PCA: Pearson, K. (1901) On lines and planes of closest fit to points in space. Philosophical Magazine, 2, 559–572.
#' * PLS: Wold, S., Ruhe, A., Wold, H., and Dunn, W.J. (1984) The collinearity problem in linear regression. the partial least squares (PLS) approach to generalized inverses. SIAM Journal of Scientific and Statistical Computing, 5, 735–743.
#' * CCA: Hotelling, H. (1936) Relations between two sets of variates. Biometrika, 28, 321–377.
#' * IFA: Tucker, L. R. (1958). An inter-battery method of factor analysis. Psychometrika, 23(2), 111-136.
#' * GSVD: Van Loan, C. (1976) Generalizing the singular value decomposition. SIAM Journal on Numerical Analysis, 13, 76–83.
#' 
#' @examples 
#' data(potato)
#' X <- potato$Chemical
#' y <- potato$Sensory[,1,drop=FALSE]
#' 
#' pca.pot  <- pca(X, ncomp = 2)
#' pcr.pot  <- pcr(y ~ X, ncomp = 2)
#' pls.pot  <- plsr(y ~ X, ncomp = 2)
#' cca.pot  <- cca(potato[1:2])
#' ifa.pot  <- ifa(potato[1:2])
#' gsvd.pot <- gsvd(lapply(potato[3:4], t))
#' 
#' @seealso Overviews of available methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' NULL

#' Principal Component Analysis - PCA
#' @param X \code{matrix} of input data.
#' @param scale \code{logical} indicating if variables should be standardised (default=FALSE).
#' @param ncomp \code{integer} number of principal components to return.
#' @param ... additional arguments to \code{pls:pcr}.
#'
#' @return \code{multiblock} object with scores, loadings, mean X values and explained variances.
#' @references Pearson, K. (1901) On lines and planes of closest fit to points in space. Philosophical Magazine, 2, 559–572.
#' @seealso Overviews of available methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' @export
pca <- function(X, scale=FALSE, ncomp=1, ...){
  X <- as.matrix(unclass(X))
  if(!inherits(X,'matrix'))
    stop("'X' must be a matrix")
  dat <- list(y=rnorm(nrow(X)), X = I(X))
  PCR <- pls::pcr(y ~ X, ncomp = ncomp, data = dat, scale = scale, ...)
  mod <- list(loadings=PCR$loadings, scores=PCR$scores, Xmeans=PCR$Xmeans, explVar <- PCR$Xvar/PCR$Xtotvar, PCA = PCR)
  mod$info <- list(method = "Principal Component Analysis", 
                   scores = "Scores", loadings = "Loadings",
                   blockScores = "Not used", blockLoadings = "Not used")
  mod$call <- match.call()
  class(mod) <- c('multiblock','list')
  return(mod)
}

#' Canonical Correlation Analysis - CCA
#' @param X \code{matrix} of input data.
#' @return \code{multiblock} object with associated with printing, scores, loadings, ...
#' @importFrom stats cancor
#' @references Hotelling, H. (1936) Relations between two sets of variates. Biometrika, 28, 321–377.
#' @seealso Overviews of available methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' @export
cca <- function(X){
  cc <- cancor(X[[1]], X[[2]])
  loadings <- cc[c('xcoef','ycoef')]
  scores   <- list((X[[1]]-rep(cc$xcenter, each=nrow(X[[1]])))%*%loadings[[1]],
                   (X[[2]]-rep(cc$ycenter, each=nrow(X[[2]])))%*%loadings[[2]])
  names(loadings) <- names(scores) <- names(X)
  mod <- list(blockLoadings=loadings, blockScores=scores, CCA = cc)
  mod$info <- list(method = "Canonical Correlation Analysis", 
                   scores = "Not used", loadings = "Not used",
                   blockScores = "Projected blocks", blockLoadings = "Block loadings")
  mod$call <- match.call()
  class(mod) <- c('multiblock','list')
  return(mod)
}


#' Inter-battery Factor Analysis - IFA
#' @param X \code{matrix} of input data.
#' @param ncomp \code{integer} number of principal components to return.
#' @param scale \code{logical} indicating if variables should be standardised (default=FALSE).
#' @param verbose \code{logical} indicating if intermediate results should be printed.
#' @param ... additional arguments to \code{pls:pcr}.
#'
#' @return \code{multiblock} object with associated with printing, scores, loadings, ...
#' @references Tucker, L. R. (1958). An inter-battery method of factor analysis. Psychometrika, 23(2), 111-136.
#' @seealso Overviews of available methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' @export
ifa <- function(X, ncomp=1, scale=FALSE, verbose=FALSE, ...){
  # InterbatteryFactorAnalysis
  # Only for two blocks
  X <- lapply(X, function(i) scale(i, scale = FALSE))
  if(length(ncomp)==1){ ncomp <- rep(ncomp,length(X)) }
  res <- RGCCA::rgcca(A = X, C = matrix(c(0,1,1,0),2,2), tau=c(1,1), verbose = verbose, scale = scale, ncomp=ncomp, ...)
  # A = coef. left, B = coef. right, corr = canonical correlation
  loadings <- list(res$astar[[1]], res$astar[[2]])
  scores   <- list(X[[1]]%*%loadings[[1]], X[[2]]%*%loadings[[2]])
  names(scores) <- names(loadings) <- names(X)
  corr <- sqrt(res$AVE$AVE_inner)
  mod <- list(blockLoadings=loadings, blockScores=scores, corr=corr, X = X, rgcca = res)
  mod$info <- list(method = "Interbattery Factor Analysis", 
                   scores = "Not used", loadings = "Not used",
                   blockScores = "Projected blocks", blockLoadings = "Block loadings")
  mod$call <- match.call()
  class(mod) <- c('multiblock','list')
  return(mod)
}


#' Generalised Singular Value Decomposition
#' @param X \code{matrix} of input data.
#' @return \code{multiblock} object with associated with printing, scores, loadings, ...
#' @importFrom stats cancor
#' @references Van Loan, C. (1976) Generalizing the singular value decomposition. SIAM Journal on Numerical Analysis, 13, 76–83.
#' @seealso Overviews of available methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' @export
gsvd  <- function(X){
  res <- geigen::gsvd(as.matrix(X[[1]]), as.matrix(X[[2]]))
  loadings <- res$Q
  scores   <- list(res$U, res$V)
  names(scores) <- names(X)
  mod <- list(loadings=loadings, blockScores=scores, GSVD = res)
  mod$info <- list(method = "Generalized Singular Value Decomposition", 
                   scores = "Not used", loadings = "Loadings",
                   blockScores = "Block scores", blockLoadings = "Not used")
  mod$call <- match.call()
  class(mod) <- c('multiblock','list')
  return(mod)
}
