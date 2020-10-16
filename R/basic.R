#' @name basic
#' @aliases pca pcr plsr ifa
#' @title Single- and Two-Block Methods
#' @description This documentation covers a range of single- and two-block methods. In particular:
#' \itemize{
#' \item Principal Component Analysis (\code{pca})  
#' \item Principal Component Regression (\code{pcr})  
#' \item Partial Least Squares Regression (\code{plsr})  
#' \item Canonical Correlation Analysis (\code{cca})  
#' \item Interbattery Factor Analysis (\code{ifa})
#' }  
#' 
#' @importFrom pls pcr plsr
#' @rdname basic 
#' @export
pca <- function(X, scale=FALSE, ncomp=1, ...){
  X <- as.matrix(unclass(X))
  if(!inherits(X,'matrix'))
    stop("'X' must be a matrix")
  dat <- list(y=rnorm(nrow(X)), X = I(X))
  PCR <- pls::pcr(y ~ X, ncomp = ncomp, data = dat, scale = scale)
  PCA <- PCR[c('scores','loadings','Xmeans')]
  PCA$explVar <- PCR$Xvar/PCR$Xtotvar
  PCA
  # X <- scale(X, scale = FALSE)
  # if(length(ncomp)==1){ ncomp <- rep(ncomp,length(X)) }
  # res <- RGCCA::rgcca(A = list(X,X), C = matrix(c(0,1,1,0),2,2), tau=c(1,1), verbose = verbose, scale = scale, ncomp=ncomp, ...)
  # # A = scores, B = loadings
  # return(list(A = X%*%res$astar[[1]], B = res$astar[[2]], X = X, rgcca = res))
}

#' @rdname basic 
#' @export
pls::pcr

#' @rdname basic 
#' @export
pls::plsr

#' @importFrom stats cancor
#' @rdname basic 
#' @export
cca <- function(X){
  cancor(X[1:2])
  # 
  # X <- lapply(X, function(i) scale(i, scale = FALSE))
  # if(length(ncomp)==1){ ncomp <- rep(ncomp,length(X)) }
  # res <- RGCCA::rgcca(A = X, C = matrix(c(0,1,1,0),2,2), tau=c(0,0), verbose = verbose, scale = scale, ncomp=ncomp, ...)
  # # A = coef. left, B = coef. right, corr = canonical correlation
  # return(list(A = res$astar[[1]], B = res$astar[[2]], corr = sqrt(res$AVE$AVE_inner), X = X, rgcca = res))
}

#' @rdname basic 
#' @export
ifa <- function(X, ncomp=1, scale=FALSE, verbose=FALSE, ...){
  # InterbatteryFactorAnalysis
  # Only for two blocks
  X <- lapply(X, function(i) scale(i, scale = FALSE))
  if(length(ncomp)==1){ ncomp <- rep(ncomp,length(X)) }
  res <- RGCCA::rgcca(A = X, C = matrix(c(0,1,1,0),2,2), tau=c(1,1), verbose = verbose, scale = scale, ncomp=ncomp, ...)
  # A = coef. left, B = coef. right, corr = canonical correlation
  return(list(A = res$astar[[1]], B = res$astar[[2]], corr = sqrt(res$AVE$AVE_inner), X = X, rgcca = res))
}

