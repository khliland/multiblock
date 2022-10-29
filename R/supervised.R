#' @name supervised
#' @title Supervised Multiblock Methods
#' 
#' @description Collection of supervised multiblock methods:
#' * MB-PLS - Multiblock Partial Least Squares (\code{\link{mbpls}})
#' * sMB-PLS - Sparse Multiblock Partial Least Squares (\code{\link{smbpls}})
#' * SO-PLS - Sequential and Orthogonalized PLS (\code{\link{sopls}})
#' * PO-PLS - Parallel and Orthogonalized PLS (\code{\link{popls}})
#' * ROSA - Response Oriented Sequential Alternation (\code{\link{rosa}})
#' * mbRDA - Multiblock Redundancy Analysis (\code{\link{mbrda}})
#' 
#' @importFrom RGCCA rgcca
#' @importFrom ade4 mbpcaiv ktab.list.df dudi.pca
#' 
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results and plotting are found in \code{\link{multiblock_results}} and \code{\link{multiblock_plots}}, respectively.
#'  
#' @examples
#' data(potato)
#' mb <- mbpls(Sensory ~ Chemical + Compression, data=potato, ncomp = 5)
#' print(mb)
#' 
#' # Convert data.frame with AsIs objects to list of matrices
#' potatoList <- lapply(potato, unclass)
#' mbr <- mbrda(Sensory ~ Chemical + Compression, data=potatoList, ncomp = 10)
#' print(mbr)
#' scoreplot(mbr, labels="names")
#' 
NULL

#' Multiblock Partial Least Squares - MB-PLS
#' 
#' @param formula Model formula accepting a single response (block) and predictor block names separated by + signs.
#' @param data The data set to analyse.
#' @param subset Expression for subsetting the data before modelling.
#' @param na.action How to handle NAs (no action implemented).
#' @param X \code{list} of input blocks. If X is supplied, the formula interface is skipped.
#' @param Y \code{matrix} of responses.
#' @param ncomp \code{integer} number of PLS components.
#' @param scale \code{logical} for autoscaling inputs (default = FALSE).
#' @param blockScale Either a \code{character} indicating type of block scaling or a \code{numeric} vector of block weights (see Details).
#' @param ... additional arguments to pls::plsr.
#' 
#' @description A function computing MB-PLS scores, loadings, etc. on the super-level and 
#' block-level.
#' 
#' @details MB-PLS is the prototypical component based supervised multiblock method.
#' It was originally formulated as a two-level method with a block-level and a super-level,
#' but it was later discovered that it could be expressed as an ordinary PLS on concatenated
#' weighted X blocks followed by a simple loop for calculating block-level loading weights,
#' loadings and scores. This implementation uses the \code{\link[pls]{plsr}} function on the
#' scaled input blocks (1/sqrt(ncol)) enabling all summaries and plots from the \code{pls}
#' package.
#' 
#' Block weighting is performed after scaling all variables and is by default 
#' \code{"sqrtnvar"}: 1/sqrt(ncol(X\[\[i\]\])) in each block. Alternatives
#' are \code{"ssq"}: 1/norm(X\[\[i\]\], "F")^2 and \code{"none"}: 1/1. Finally, if
#' a \code{numeric} vector is supplied, it will be used to scale the blocks
#' after \code{"ssq"} scaling, i.e., Z\[\[i\]\] = X\[\[i\]\] / norm(X\[\[i\]\], "F")^2 * blockScale\[i\].
#' 
#' @return \code{multiblock, mvr} object with super-scores, super-loadings, block-scores and block-loading, and the underlying 
#' \code{mvr} (PLS) object for the super model, with all its result and plot possibilities. Relevant plotting functions: \code{\link{multiblock_plots}} 
#' and result functions: \code{\link{multiblock_results}}.
#' 
#' @references 
#' * Wangen, L.E. and Kowalski, B.R. (1988). A multiblock partial least squares algorithm for investigating complex chemical systems. Journal of Chemometrics, 3, 3–20.
#' * Westerhuis, J.A., Kourti, T., and MacGregor,J.F. (1998). Analysis of multiblock and hierarchical PCA and PLS models. Journal of Chemometrics, 12, 301–321.
#' 
#' @examples 
#' data(potato)
#' # Formula interface
#' mb <- mbpls(Sensory ~ Chemical+Compression, data=potato, ncomp = 5)
#' 
#' # ... or X and Y
#' mb.XY <- mbpls(X=potato[c('Chemical','Compression')], Y=potato[['Sensory']], ncomp = 5)
#' identical(mb$scores, mb.XY$scores)
#' print(mb)
#' scoreplot(mb, labels="names") # Exploiting mvr object structure from pls package
#' 
#' # Block scaling with emphasis on first block
#' mbs <- mbpls(Sensory ~ Chemical+Compression, data=potato, ncomp = 5, blockScale = c(10, 1))
#' scoreplot(mbs, labels="names") # Exploiting mvr object structure from pls package
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' @export
mbpls <- function(formula, data, subset, na.action, X=NULL, Y=NULL, ncomp=1, scale=FALSE, blockScale=c("sqrtnvar","ssq","none"), ...){
  if(is.null(X)){ # Use formula interface is X is not supplied.
    ## Get the model frame
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]                # Retain only the named arguments
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    X <- mf[-1]
    
    ## Get the terms
    mt <- attr(mf, "terms")        # This is to include the `predvars'
    # attribute of the terms
    
    ## Get the data matrices
    Y <- model.response(mf, "numeric")
    if (is.matrix(Y)) {
      if (is.null(colnames(Y))) 
        colnames(Y) <- paste("Y", 1:dim(Y)[2], sep = "")
    }
    else {
      Y <- as.matrix(Y)
      colnames(Y) <- deparse(formula[[2]])
    }
  }
  
  # Block scaling
  nblock <- length(X)
  Xo <- lapply(X, function(x)if(is.factor(x)){return(dummycode(x))}else{return(x)})
  Xo <- lapply(lapply(Xo, as.matrix), function(x)scale(x, scale=scale))
  if(is.numeric(blockScale)){ # User supplied scale
    X <- lapply(1:nblock, function(i) Xo[[i]]/norm(Xo[[i]],"F")^2*blockScale[i])
  } else {
    if(is.character(blockScale)){
      if(blockScale[1] == "sqrtnvar")
        X <- lapply(Xo, function(x) x/sqrt(ncol(x)))
      else if(blockScale[1] == "ssq")
        X <- lapply(1:nblock, function(i) Xo[[i]]/norm(Xo[[i]],"F")^2)
      else if(blockScale[1] == "none")
        X <- Xc
      else
        stop("Unknown format for 'blockScale'")
    } else {
      stop("Unknown format for 'blockScale'")
    }
  }
  Y  <- as.matrix(Y)
  Xc  <- do.call(cbind, X)
  dat <- list(X = Xc, Y = Y)
  comps <- paste('Comp', 1:ncomp)
  mod <- pls::plsr(Y ~ X, data = dat, ncomp = ncomp, ...)
  U   <- normCols(mod$Yscores)^2 # normalized Y-scores
  Wb  <- Tb <- Pb <- list()
  Wt  <- matrix(0, length(X), ncomp)
  for(b in 1:length(X)){ # Loop over blocks
    Wb[[b]] <- crossprod(Xo[[b]], U)             # Block loading weights
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
  mod$blockScale <- blockScale
  mod$data <- list(X = X, Y = Y)
  attr(mod$scores, "explvar") <- attr(mod$loadings, "explvar") <- mod$Xvar/mod$Xtotvar*100
  mod$explvar <- mod$Xvar/mod$Xtotvar*100
  # for(i in 1:length(X)){
  #   expli <- numeric(ncomp)
  #   normi <- base::norm(X[[i]],'F')^2
  #   for(j in 1:ncomp){
  #     expli[j] <- 1-base::norm(X[[i]]-mod$blockScores[[i]][,1:j,drop=FALSE] %*% t(mod$blockLoadings[[i]][,1:j,drop=FALSE]), 'F')^2/normi*100
  #   }
  #   attr(mod$blockScores[[i]], "explvar") <- attr(mod$blockLoadings[[i]], "explvar") <- diff(c(0,expli))
  # }
  mod$info <- list(method = "Multiblock PLS", 
                   scores = "Superscores", loadings = "Superloadings",
                   blockScores = "Block scores", blockLoadings = "Block loadings")
  mod$call <- match.call()
  class(mod) <- c('multiblock','mvr')
  return(mod)
}

#' Sparse Multiblock Partial Least Squares - sMB-PLS
#' 
#' @param formula Model formula accepting a single response (block) and predictor block names separated by + signs.
#' @param data The data set to analyse.
#' @param subset Expression for subsetting the data before modelling.
#' @param na.action How to handle NAs (no action implemented).
#' @param X \code{list} of input blocks. If X is supplied, the formula interface is skipped.
#' @param Y \code{matrix} of responses.
#' @param ncomp \code{integer} number of PLS components.
#' @param scale \code{logical} for autoscaling inputs (default = FALSE).
#' @param shrink \code{numeric} scalar indicating degree of L1-shrinkage/Soft-Thresholding (optional), 0 <= shrink < 1.
#' @param truncation \code{character} indicating type of truncation (optional) "Lenth" uses 
#' asymmetric confidence intervals to determine outlying loading weights. "quantile" uses
#' a quantile plot approach to determining outliers.
#' @param trunc.width \code{numeric} indicating confidence of "Lenth type" confidence interval
#' or quantile in "quantile plot" approach. Default = 0.95.
#' @param blockScale Either a \code{character} indicating type of block scaling or a \code{numeric} vector of block weights (see Details).
#' @param ... additional arguments to pls::plsr.
#' 
#' @description sMB-PLS is an adaptation of MB-PLS (\code{\link{mbpls}}) that enforces sparseness in loading weights
#' when computing PLS components in the global model. 
#' 
#' @details Two versions of sparseness are supplied: Soft-Threshold PLS, also
#' known as Sparse PLS, and Truncation PLS. The former uses L1 shrinkage of loading weights, while the latter
#' comes in two flavours, both estimating inliers and outliers. The "Lenth" method uses asymmetric confidence
#' intervals around the median of a loading weigh vector to estimate inliers. The "quantile" method uses
#' a quantile plot approach to estimate outliers as deviations from the estimated quantile line. As with 
#' ordinary MB-PLS scaled input blocks (1/sqrt(ncol)) are used.
#' 
#' Block weighting is performed after scaling all variables and is by default 
#' \code{"sqrtnvar"}: 1/sqrt(ncol(X\[\[i\]\])) in each block. Alternatives
#' are \code{"ssq"}: 1/norm(X\[\[i\]\], "F")^2 and \code{"none"}: 1/1. Finally, if
#' a \code{numeric} vector is supplied, it will be used to scale the blocks
#' after \code{"ssq"} scaling, i.e., Z\[\[i\]\] = X\[\[i\]\] / norm(X\[\[i\]\], "F")^2 * blockScale\[i\].
#' 
#' @return \code{multiblock, mvr} object with super-scores, super-loadings, block-scores and block-loading, and the underlying 
#' \code{mvr} (PLS) object for the super model, with all its result and plot possibilities. Relevant plotting functions: \code{\link{multiblock_plots}} 
#' and result functions: \code{\link{multiblock_results}}.
#' 
#' @references 
#' * Sæbø, S.; Almøy, T.; Aarøe, J. & Aastveit, A. ST-PLS: a multi-directional nearest shrunken centroid type classifier via PLS Journal of Chemometrics: A Journal of the Chemometrics Society, Wiley Online Library, 2008, 22, 54-62.
#' * Lê Cao, K.; Rossouw, D.; Robert-Granié, C. & Besse, P. A sparse PLS for variable selection when integrating omics data Statistical applications in genetics and molecular biology, 2008, 7.
#' * Liland, K.; Høy, M.; Martens, H. & Sæbø, S. Distribution based truncation for variable selection in subspace methods for multivariate regression Chemometrics and Intelligent Laboratory Systems, 2013, 122, 103-111.
#' * Karaman, I.; Nørskov, N.; Yde, C.; Hedemann, M.; Knudsen, K. & Kohler, A. Sparse multi-block PLSR for biomarker discovery when integrating data from LC--MS and NMR metabolomics Metabolomics, 2015, 11, 367-379.
#' 
#' @examples 
#' data(potato)
#' 
#' # Truncation MB-PLS 
#' # Loading weights inside 60% confidence intervals around the median are set to 0.
#' tmb <- smbpls(Sensory ~ Chemical+Compression, data=potato, ncomp = 5, 
#'               truncation = "Lenth", trunc.width = 0.6)
#'               
#' # Alternative XY-interface
#' tmb.XY <- smbpls(X=potato[c('Chemical','Compression')], Y=potato[['Sensory']], ncomp = 5, 
#'               truncation = "Lenth", trunc.width = 0.6)
#' identical(tmb, tmb.XY)
#' scoreplot(tmb, labels="names") # Exploiting mvr object structure from pls package
#' loadingweightplot(tmb, labels="names")
#' 
#' # Soft-Threshold / Sparse MB-PLS 
#' # Loading weights are subtracted by 60% of maximum value.
#' smb <- smbpls(X=potato[c('Chemical','Compression')], Y=potato[['Sensory']], 
#'               ncomp = 5, shrink = 0.6)
#' print(smb)
#' scoreplot(smb, labels="names") # Exploiting mvr object structure from pls package
#' loadingweightplot(smb, labels="names")
#' 
#' # Emphasis may be different for blocks
#' smb <- smbpls(X=potato[c('Chemical','Compression')], Y=potato[['Sensory']], 
#'               ncomp = 5, shrink = 0.6, blockScale = c(1, 10))
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' @importFrom plsVarSel truncation stpls mvrV
#' @export
smbpls <- function(formula, data, subset, na.action, X=NULL, Y=NULL, ncomp=1, scale=FALSE, shrink=NULL, truncation=NULL, trunc.width=0.95, blockScale=c("sqrtnvar","ssq","none"),...){
  if(is.null(truncation) && is.null(shrink)){
    stop("Specify either 'truncation' or 'shrink' to define sparseness method.")
  }
  if(!is.null(truncation) && !is.null(shrink)){
    stop("Choose either 'truncation' or 'shrink', not both.")
  }
  
  if(is.null(X)){ # Use formula interface is X is not supplied.
    ## Get the model frame
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]                # Retain only the named arguments
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    X <- mf[-1]

    ## Get the terms
    mt <- attr(mf, "terms")        # This is to include the `predvars'
    # attribute of the terms
    
    ## Get the data matrices
    Y <- model.response(mf, "numeric")
    if (is.matrix(Y)) {
      if (is.null(colnames(Y))) 
        colnames(Y) <- paste("Y", 1:dim(Y)[2], sep = "")
    }
    else {
      Y <- as.matrix(Y)
      colnames(Y) <- deparse(formula[[2]])
    }
  }

  # Block scaling
  nblock <- length(X)
  Xo <- lapply(X, function(x)if(is.factor(x)){return(dummycode(x))}else{return(x)})
  Xo <- lapply(lapply(Xo, as.matrix), function(x)scale(x, scale=scale))
  if(is.numeric(blockScale)){ # User supplied scale
    X <- lapply(1:nblock, function(i) Xo[[i]]/norm(Xo[[i]],"F")^2*blockScale[i])
  } else {
    if(is.character(blockScale)){
      if(blockScale[1] == "sqrtnvar")
        X <- lapply(Xo, function(x) x/sqrt(ncol(x)))
      else if(blockScale[1] == "ssq")
        X <- lapply(1:nblock, function(i) Xo[[i]]/norm(Xo[[i]],"F")^2)
      else if(blockScale[1] == "none")
        X <- Xc
      else
        stop("Unknown format for 'blockScale'")
    } else {
      stop("Unknown format for 'blockScale'")
    }
  }
  Y  <- as.matrix(Y)
  Xc  <- do.call(cbind, X)
  dat <- list(X = Xc, Y = Y)
  comps <- paste('Comp', 1:ncomp)
  if(!is.null(truncation)){
    methodName <- "Sparse Multiblock PLS (Truncation)"
    mod <- plsVarSel::truncation(Y ~ X, data = dat, ncomp = ncomp, truncation = truncation, trunc.width = trunc.width, ...)
  }
  if(!is.null(shrink)){
    if(length(shrink)>1)
      warning("Using only first value of 'shrink'")
    methodName <- "Sparse Multiblock PLS (Soft-Threshold)"
    mod <- plsVarSel::stpls(Y ~ X, data = dat, ncomp = ncomp, shrink = shrink[1], method = "stpls", ...)
    mod$Yscores         <- mod$Yscores[,,1]
    mod$scores          <- mod$scores[,,1]
    mod$loadings        <- mod$loadings[,,1]
    mod$Yloadings       <- mod$Yloadings[,,1]
    mod$loading.weights <- mod$loading.weights[,,1]
    mod$projection      <- mod$projection[,,1]
  }
  U   <- normCols(mod$Yscores)^2 # normalized Y-scores
  Wb  <- Tb <- Pb <- list()
  Wt  <- matrix(0, length(X), ncomp)
  for(b in 1:length(X)){ # Loop over blocks
    Wb[[b]] <- crossprod(Xo[[b]], U)             # Block loading weights
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
  mod$data <- list(X = X, Y = Y)
  attr(mod$scores, "explvar") <- attr(mod$loadings, "explvar") <- attr(mod$loading.weights, "explvar") <- mod$Xvar/mod$Xtotvar*100
  mod$explvar <- mod$Xvar/mod$Xtotvar*100

  mod$info <- list(method = methodName, 
                   scores = "Superscores", loadings = "Superloadings",
                   blockScores = "Block scores", blockLoadings = "Block loadings")
  mod$call <- match.call()
  class(mod) <- c('multiblock','mvr')
  return(mod)
}

#' Multiblock Redundancy Analysis - mbRDA
#' @param formula Model formula accepting a single response (block) and predictor block names separated by + signs.
#' @param data The data set to analyse.
#' @param subset Expression for subsetting the data before modelling.
#' @param na.action How to handle NAs (no action implemented).
#' @param X \code{list} of input blocks.
#' @param Y \code{matrix} of responses.
#' @param ncomp \code{integer} number of PLS components.
#' @param ... additional arguments to ade4::mbpcaiv.
#' 
#' @description This is a wrapper for the \code{ade4::mbpcaiv} function for computing mbRDA.
#' 
#' @details mbRDA is a multiblock formulation of Redundancy (Data) Analysis. RDA is theoretically
#' between PLS and GCA. Like GCA, RDA does not consider correlations within X, but like
#' PLS it does consider correlations within Y. RDA can also be viewed as a PCR of Y constrained to
#' have scores that are also linear combinations of X. If the \code{adegraphics} package is attached,
#' a nice overview can be plotted as \code{plot(mbr$mbpcaiv)} following the example below.
#' 
#' @references Bougeard, S., Qannari, E.M., Lupo, C., andHanafi, M. (2011). From Multiblock Partial Least Squares to Multiblock Redundancy Analysis. A Continuum Approach. Informatica, 22(1), 11–26.
#' 
#' @return \code{multiblock, mvr} object with scores, block-scores and block-loading. Relevant plotting functions: \code{\link{multiblock_plots}} 
#' and result functions: \code{\link{multiblock_results}}.
#' 
#' @examples
#' # Convert data.frame with AsIs objects to list of matrices
#' data(potato)
#' potatoList <- lapply(potato, unclass)
#' 
#' mbr <- mbrda(Sensory ~ Chemical + Compression, data = potatoList, ncomp = 10)
#' mbr.XY <- mbrda(X = potatoList[c('Chemical','Compression')], Y = potatoList[['Sensory']], 
#'                 ncomp = 10)
#' print(mbr)
#' scoreplot(mbr) # Exploiting mvr object structure from pls package
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' @export
mbrda <- function(formula, data, subset, na.action, X=NULL, Y=NULL, ncomp=1, ...){
  # MBRedundancyAnalysis
  if(is.null(X)){ # Use formula interface is X is not supplied.
    ## Get the model frame
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]                # Retain only the named arguments
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    X <- mf[-1]
    
    ## Get the terms
    mt <- attr(mf, "terms")        # This is to include the `predvars'
    # attribute of the terms
    
    ## Get the data matrices
    Y <- model.response(mf, "numeric")
    if (is.matrix(Y)) {
      if (is.null(colnames(Y))) 
        colnames(Y) <- paste("Y", 1:dim(Y)[2], sep = "")
    }
    else {
      Y <- as.matrix(Y)
      colnames(Y) <- deparse(formula[[2]])
    }
  }
  
  dat <- list(X = X, Y = Y)
  Xd  <- lapply(X,as.data.frame)
  Xk  <- ktab.list.df(Xd)
  Y   <- dudi.pca(as.data.frame(Y), scannf = FALSE, nf = ncol(Y))
  res <- mbpcaiv(Y, Xk, scale = TRUE, scannf = FALSE, nf = ncomp, ...)
  
  varT <- diag(crossprod(res$lX * res$lw, res$lX))
  covarTY <- diag(tcrossprod(crossprod(res$lX * res$lw, 
                                       as.matrix(res$tabY))))
  varExplTY <- (covarTY/varT)/sum(covarTY/varT) * 100
  vars <- cumsum(unlist(lapply(X, ncol))); v <- 1
  blockScores <- blockLoadings <- list()
  for(i in 1:length(vars)){
    blockLoadings[[i]] <- res$Tfa[v:vars[i],]
    blockScores[[i]] <- res$Tl1[(i-1)*nrow(X[[i]])+(1:nrow(X[[i]])),]
    v <- vars[i]+1
  }
  names(blockLoadings) <- names(X)
  for(i in 1:length(X)){
    rownames(blockLoadings[[i]]) <- colnames(X[[i]])
    rownames(blockScores[[i]])   <- rownames(X[[i]])
    colnames(blockScores[[i]])   <- colnames(blockLoadings[[i]]) <- paste0("Comp ", 1:ncomp)
  }
  #                     u                 v
  mod <- list(Yscores=res$lY, Yloadings=res$Yc1, scores=res$lX, 
              blockLoadings=blockLoadings, blockScores=blockScores, 
              varT=varT, covarTY=covarTY, varExplTY=varExplTY, data = dat, mbpcaiv=res)
  mod$info <- list(method = "Multiblock RDA", 
                   scores = "Scores", loadings = "Not used",
                   blockScores = "Block scores", blockLoadings = "Block loadings") # t_m, w_m
  mod$call <- match.call()
  class(mod) <- c('multiblock','mvr')
  return(mod)
}

