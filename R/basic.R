#' @name basic
#' @aliases pca pcr plsr cca gsvd ifa
#' @title Single- and Two-Block Methods
#' @description This documentation covers a range of single- and two-block methods. In particular:
#' * PCA - Principal Component Analysis (\code{pca})  
#' * PCR - Principal Component Regression (\code{pcr})  
#' * PLSR - Partial Least Squares Regression (\code{plsr})  
#' * CCA - Canonical Correlation Analysis (\code{cca})  
#' * IFA - Interbattery Factor Analysis (\code{ifa})
#' * GSVD - Generalized SVD (\code{gsvd})
#' }  
#' 
#' @importFrom pls pcr plsr
#' @importFrom geigen gsvd
#' @rdname basic 
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

#' @rdname basic 
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


#' @rdname basic 
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
