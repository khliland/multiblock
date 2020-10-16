#' @name supervised
#' @title Supervised Multiblock Methods
#' @aliases mbpls mbra
#' 
#' @description Collection of supervised multiblock methods:
#' * MB-PLS - Multiblock Partial Least Squares (\code{mbpls})
#' * SO-PLS - Sequential and Orthgonalized PLS (\code{\link{sopls}})
#' * ROSA - Response Oriented Sequential Alternation (\code{\link{rosa}})
#' * Multiblock Redundancy Analysis (\code{mbra})
#' 
#' @importFrom RGCCA rgcca
#' @importFrom ade4 mbpcaiv ktab.list.df dudi.pca
# #' @importFrom MetStaT ASCA.Calculate
#'  

#' @rdname supervised 
#' @export
mbpls <- function(X, Y, ncomp=1, scale=FALSE, ...){
  # TODO: Extend with various block norms
  X <- lapply(lapply(X, as.matrix), function(x) x/sqrt(ncol(x)))
  Y <- as.matrix(Y)
  Xc <- do.call(cbind, X)
  dat <- list(X = Xc, Y = Y)
  mod <- pls::plsr(Y ~ X, data=dat, ncomp=ncomp, scale=scale, ...)
  U <- normCols(mod$Yscores)^2 # normalized Y-scores
  Wb <- Tb <- Pb <- list()
  Wt <- matrix(0, length(X), ncomp)
  for(b in 1:length(X)){ # Loop over blocks
    Wb[[b]] <- crossprod(X[[b]], U)                    # Block weights
    Tb[[b]] <- X[[b]] %*% Wb[[b]] / sqrt(ncol(X[[b]])) # Block scores
    Pb[[b]] <- crossprod(X[[b]], Tb[[b]])              # Block loadings
    for(k in 1:ncomp){
      Pb[[b]][,k] <- Pb[[b]][,k]/drop(crossprod(Tb[[b]][,k]))
    }
    Wt[b,] <- colSums(crossprod(Tb[[b]], U)) # Super weights
  }
  names(Tb) <- names(Pb) <- names(Wb) <- names(X)
  return(list(Tt=mod$scores, Wt=Wt, Tb=Tb, Wb=Wb, Pb=Pb, plsmod=mod))
}

#' @rdname supervised 
#' @export
mbra <- function(X, Y, ncomp=1, ...){
  # MBRedundancyAnalysis
  Xd <- lapply(X,as.data.frame)
  Xk <- ktab.list.df(Xd)
  Y  <- dudi.pca(Y, scannf = FALSE, nf = 1)
  res <- mbpcaiv(Y, Xk, scale = TRUE, scannf = FALSE, nf = ncomp, ...)
  
  varT <- diag(crossprod(res$lX * res$lw, res$lX))
  covarTY <- diag(tcrossprod(crossprod(res$lX * res$lw, 
                                       as.matrix(res$tabY))))
  varExplTY <- (covarTY/varT)/sum(covarTY/varT) * 100
  return(list(Yscores=res$lY, Yloadings=res$Yc1, scores=res$lX, loadings=res$Tfa, varT=varT, covarTY=covarTY, varExplTY=varExplTY, mbpcaivObject=res))
}

# #' @rdname supervised 
# #' @export
# ASCA <- function(X, ...){
#   # ASCA.Calculate(data, levels, equation.elements = "", scaling = FALSE, 
#   #                only.means.matrix = FALSE, use.previous.asca = NULL)
#   "MetStaT (not model ellipsoids)"}

