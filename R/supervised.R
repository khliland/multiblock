#' @name supervised
#' @title Supervised Multiblock Methods
#' 
#' @description Collection of supervised multiblock methods:
#' * MB-PLS - Multiblock Partial Least Squares (\code{\link{mbpls}})
#' * SO-PLS - Sequential and Orthogonalized PLS (\code{\link{sopls}})
#' * PO-PLS - Parallel and Orthogonalized PLS (\code{\link{popls}})
#' * ROSA - Response Oriented Sequential Alternation (\code{\link{rosa}})
#' * mbRDA - Multiblock Redundancy Analysis (\code{\link{mbrda}})
#' 
#' @importFrom RGCCA rgcca
#' @importFrom ade4 mbpcaiv ktab.list.df dudi.pca
#' 
#' @seealso Overviews of available methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
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
#' scoreplot(mbr, labels="names")
#' 
NULL

#' Multiblock Partial Least Squares - MB-PLS
#' 
#' @param X \code{list} of input blocks.
#' @param Y \code{matrix} of responses.
#' @param ncomp \code{integer} number of PLS components.
#' @param scale \code{logical} for autoscaling inputs (default = FALSE).
#' @param ... additional arguments to pls::plsr.
#' 
#' @description MB-PLS is the prototypical component based supervised multiblock method.
#' It was originally formulated as a two-level method with a block-level and a super-level,
#' but it was later discovered that it could be expressed as an ordinary PLS on concatenated
#' weighted X blocks followed by a simple loop for calculating block-level loading weights,
#' loadings and scores. This implementation uses the \code{\link[pls]{plsr}} function on the
#' scaled input blocks (1/sqrt(ncol)) enabling all summaries and plots from the \code{pls}
#' package.
#' 
#' @return \code{mbpls} object containing the underlying \code{pls} object, with all its result and plot possibilities plus block-wise loadings, loading weights and scores.
#' 
#' @references 
#' * Wangen, L.E. and Kowalski, B.R. (1988). A multiblock partial least squares algorithm for investigating complex chemical systems. Journal of Chemometrics, 3, 3–20.
#' * Westerhuis, J.A., Kourti, T., and MacGregor,J.F. (1998). Analysis of multiblock and hierarchical PCA and PLS models. Journal of Chemometrics, 12, 301–321.
#' 
#' @examples 
#' data(potato)
#' mb <- mbpls(potato[c('Chemical','Compression')], potato[['Sensory']], ncomp = 5)
#' print(mb)
#' scoreplot(mb, labels="names") # Exploiting mvr object structure from pls package
#' @seealso Overviews of available methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
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
  mod$info <- list(method = "Multiblock PLS", 
                   scores = "Superscores", loadings = "Superloadings",
                   blockScores = "Block scores", blockLoadings = "Block loadings")
  mod$call <- match.call()
  class(mod) <- c('multiblock','mvr')
  return(mod)
}

#' Multiblock Redundancy Analysis - mbRDA
#' @param X \code{list} of input blocks.
#' @param Y \code{matrix} of responses.
#' @param ncomp \code{integer} number of PLS components.
#' @param ... additional arguments to ade4::mbpcaiv.
#' 
#' @return \code{mbrda,mvr} object containing elements corresponding to a \code{pls} object, with all its result and plot possibilities plus block-wise loadings, loading weights and scores.
#' 
#' @description mbRDA is a multiblock formulation of Redundancy (Data) Analysis. RDA is theoretically
#' between PLS and GCA. Like GCA, RDA does not consider correlations within X, but like
#' PLS it does consider correlations within Y. RDA can also be viewed as a PCR of Y constrained to
#' have scores that are also linear combinations of X.
#' 
#' @references Bougeard, S., Qannari, E.M., Lupo, C., andHanafi, M. (2011). From Multiblock Partial Least Squares to Multiblock Redundancy Analysis. A Continuum Approach. Informatica, 22(1), 11–26.
#' 
#' @examples
#' # Convert data.frame with AsIs objects to list of matrices
#' data(potato)
#' potatoList <- lapply(potato, unclass)
#' 
#' mbr <- mbrda(potatoList[c('Chemical','Compression')], potatoList[['Sensory']], ncomp = 10)
#' print(mbr)
#' scoreplot(mbr) # Exploiting mvr object structure from pls package
#' @seealso Overviews of available methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
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
  mod$info <- list(method = "Multiblock RDA", 
                   scores = "Scores", loadings = "Loadings",
                   blockScores = "Not used", blockLoadings = "Not used")
  mod$call <- match.call()
  class(mod) <- c('multiblock','mvr')
  return(mod)
}

