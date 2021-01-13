#' @name supervised
#' @title Supervised Multiblock Methods
#' @aliases mbpls mbrda
#' 
#' @description Collection of supervised multiblock methods:
#' * MB-PLS - Multiblock Partial Least Squares (\code{mbpls})
#' * SO-PLS - Sequential and Orthogonalized PLS (\code{\link{sopls}})
#' * ROSA - Response Oriented Sequential Alternation (\code{\link{rosa}})
#' * mbRDA - Multiblock Redundancy Analysis (\code{mbrda})
#' 
#' @importFrom RGCCA rgcca
#' @importFrom ade4 mbpcaiv ktab.list.df dudi.pca
#'  
#' @examples
#' data(potato)
#' mb <- mbpls(potato[c('Chemical','Compression')], potato[['Sensory']], ncomp = 5)
#' print(mb)
#' 
#' # Convert data.frame with AsIs objects to list of matrices
#' potatoList <- lapply(potato, unclass)
#' mbr <- mbrda(potatoList[c('Chemical','Compression')], potatoList[['Sensory']], ncomp = 10)
#' print(mbr)
#' scoreplot(mbr)
#' 
#' TODO: 
#' * mbRDA to formula interface
#' * MBPLS to formula interface
#' * sparse MBPLS (ST-type, Lenth type)
#' * PO-PLS
#' 
#' @export
mbpls <- function(X, Y, ncomp=1, scale=FALSE, ...){
  # TODO: Extend with various block norms
  X <- lapply(lapply(X, as.matrix), function(x)scale(x, scale=scale))
  X <- lapply(X, function(x) x/sqrt(ncol(x)))
  Y <- as.matrix(Y)
  Xc  <- do.call(cbind, X)
  dat <- list(X = Xc, Y = Y)
  comps <- paste('Comp', 1:ncomp)
  mod <- pls::plsr(Y ~ X, data=dat, ncomp=ncomp, ...)
  U   <- normCols(mod$Yscores)^2 # normalized Y-scores
  Wb  <- Tb <- Pb <- list()
  Wt  <- matrix(0, length(X), ncomp)
  for(b in 1:length(X)){ # Loop over blocks
    Wb[[b]] <- crossprod(X[[b]], U)             # Block loading weights
    for(k in 1:ncomp){
      Wb[[b]][,k] <- Wb[[b]][,k]/sqrt(drop(crossprod(Wb[[b]][,k])))
    }
    Tb[[b]] <- X[[b]] %*% Wb[[b]]               # Block scores
    Pb[[b]] <- crossprod(X[[b]], Tb[[b]])       # Block loadings
    for(k in 1:ncomp){
      Pb[[b]][,k] <- Pb[[b]][,k]/drop(crossprod(Tb[[b]][,k]))
    }
    Wt[b,] <- colSums(crossprod(Tb[[b]], U)) # Super weights
    dimnames(Pb[[b]]) <- dimnames(Wb[[b]]) <- list(colnames(X[[b]]), comps)
    dimnames(Tb[[b]]) <- list(rownames(X[[b]]), comps)
  }
  dimnames(Wt) <- list(names(X), comps)
  names(Tb) <- names(Pb) <- names(Wb) <- names(X)
  mod$blockScores   <- Tb
  mod$blockLoadingweights <- Wb
  mod$blockLoadings <- Pb
  mod$superWeights  <- Wt
  mod$call <- match.call()
  class(mod) <- c('mbpls','mvr')
  # obj <- list(Tt=mod$scores, Wt=Wt, Tb=Tb, Wb=Wb, Pb=Pb, plsmod=mod, call=match.call())
  # class(obj) <- 'mbpls'
  return(mod)
}

#' @rdname supervised 
#' @export
mbrda <- function(X, Y, ncomp=1, ...){
  # MBRedundancyAnalysis
  Xd  <- lapply(X,as.data.frame)
  Xk  <- ktab.list.df(Xd)
  Y   <- dudi.pca(Y, scannf = FALSE, nf = 1)
  res <- mbpcaiv(Y, Xk, scale = TRUE, scannf = FALSE, nf = ncomp, ...)
  
  varT <- diag(crossprod(res$lX * res$lw, res$lX))
  covarTY <- diag(tcrossprod(crossprod(res$lX * res$lw, 
                                       as.matrix(res$tabY))))
  varExplTY <- (covarTY/varT)/sum(covarTY/varT) * 100
  mod <- list(Yscores=res$lY, Yloadings=res$Yc1, scores=res$lX, loadings=res$Tfa, varT=varT, covarTY=covarTY, varExplTY=varExplTY, mbpcaivObject=res)
  mod$call <- match.call()
  class(mod) <- c('mbrda','mvr')
  return(mod)
}

