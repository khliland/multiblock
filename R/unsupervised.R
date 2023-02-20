# #' @aliases sca gca mfa pcagca disco jive statis hogsvd hpca mcoa
#' @name unsupervised
#' @title Unsupervised Multiblock Methods
#' @importFrom RGCCA rgcca
#' @importFrom FactoMineR MFA GPA
#' @importFrom RegularizedSCA DISCOsca
#' @importFrom ade4 statis ktab.within withinpca
#' @importFrom r.jive jive summary.jive plot.jive
#' 
#' @description Collection of unsupervised multiblock methods:
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
#' @details 
#' Original documentation of STATIS: \link[ade4]{statis}.
#' JIVE, STATIS and HOGSVD assume variable linked matrices/data.frames, while SCA handles both links.
#' 
#' @examples 
#' # Object linked data
#' data(potato)
#' potList <- as.list(potato[c(1,2,9)])
#' pot.sca    <- sca(potList)
#' 
#' # Variable linked data
#' data(candies)
#' candyList <- lapply(1:nlevels(candies$candy),function(x)candies$assessment[candies$candy==x,])
#' can.statis <- statis(candyList)
#' plot(can.statis$statis)
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results and plotting are found in \code{\link{multiblock_results}} and \code{\link{multiblock_plots}}, respectively.
#
# Include also Penalized Exponential SCA (SLIDE, penalty-based)?????????
# 
NULL

#' Simultaneous Component Analysis - SCA
#' @param X \code{list} of input blocks.
#' @param ncomp \code{integer} number of components to extract.
#' @param scale \code{logical} indicating autoscaling of features (default = FALSE).
#' @param samplelinked \code{character/logical} indicating if blocks are linked by samples (TRUE) or variables (FALSE). Using 'auto' (default), this will be determined automatically.
#' @param ... additional arguments (not used).
#' 
#' @description This is a basic implementation of the SCA-P algorithm (least restricted SCA) with support for both
#' sample- and variable-linked modes.
#' 
#' @details SCA, in its original variable-linked version, calculates common loadings and block-wise
#' scores. There are many possible constraints and specialisations. This implementations uses
#' PCA as the backbone, thus resulting in deterministic, ordered components. A parameter controls
#' the linking mode, but if left untouched an attempt is made at automatically determining
#' variable or sample linking.
#' 
#' @return \code{multiblock} object including relevant scores and loadings. Relevant plotting functions: \code{\link{multiblock_plots}} 
#' and result functions: \code{\link{multiblock_results}}.
#' 
#' @references Levin, J. (1966) Simultaneous factor analysis of several gramian matrices. Psychometrika, 31(3), 413–419.
#' 
#' @examples
#' # Object linked data
#' data(potato)
#' potList <- as.list(potato[c(1,2,9)])
#' pot.sca    <- sca(potList)
#' plot(scores(pot.sca), labels="names")
#' 
#' # Variable linked data
#' data(candies)
#' candyList <- lapply(1:nlevels(candies$candy),function(x)candies$assessment[candies$candy==x,])
#' pot.sca    <- sca(candyList, samplelinked = FALSE)
#' pot.sca
#'
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results and plotting are found in \code{\link{multiblock_results}} and \code{\link{multiblock_plots}}, respectively.
#' @export
sca <- function(X, ncomp=2, scale=FALSE, samplelinked = 'auto', ...){
  # SVD/PCA based SVD-P with block-centring (and global scaling)
  # Infer link-mode:
  row <- unlist(lapply(X,nrow))
  col <- unlist(lapply(X,ncol))
  if(is.character(samplelinked) && samplelinked == 'auto'){
    if(all(row==row[1])){
      samplelinked <- TRUE
      if(all(col==col[1])){
        stop('Cannot determine if matrices are sample linked or not (both satisfied).')
      }
    } else {
      if(all(col==col[1])){
        samplelinked <- FALSE
      } else {
        stop('Neither rows, nor collumns match (not equal numbers in all matrices).')
      }
    }
  }
  
  # Center (and scale)
  X <- lapply(X, function(x)scale(x, scale=FALSE))
  if(samplelinked == FALSE){
    Xr <- do.call(rbind, X)
  } else {
    Xr <- do.call(cbind, X)
  }
  if(scale)
    Xr <- scale(Xr)
  # SVD/PCA
  PCA <- pca(Xr, scale=scale, ncomp=ncomp)
  if(samplelinked == FALSE){
    loadings <- PCA$loadings
    nvar   <- lapply(X, nrow)
    lookup <- unlist(lapply(1:length(X), function(i)rep(i, nvar[[i]])))
    scores <- lapply(1:length(X), function(i)PCA$scores[lookup==i,,drop=FALSE])
    names(scores) <- names(X)
    mod <- list(loadings=loadings, blockScores=scores, samplelinked=samplelinked)
    mod$info <- list(method = "Simultaneous Component Analysis",
                     scores = "Not used", loadings = "Common loadings",
                     blockScores = "Block scores", blockLoadings = "Not used")
    attr(mod$loadings, 'explvar') <- PCA$explvar
    for(i in 1:length(X)){
      explI <- numeric(length(PCA$explvar))
      denoI <- base::norm(X[[i]],'F')^2
      for(j in 1:length(PCA$explvar)){
        explI[j] <- sum(diag(loadings[,1:j,drop=FALSE] %*% tcrossprod(crossprod(scores[[i]][,1:j,drop=FALSE]), loadings[,1:j,drop=FALSE])))/denoI*100
      }
      attr(mod$blockScores[[i]],'explvar') <- diff(c(0,explI))
    }
  } else {
    scores <- PCA$scores
    nvar   <- lapply(X, ncol)
    lookup <- unlist(lapply(1:length(X), function(i)rep(i, nvar[[i]])))
    loadings <- lapply(1:length(X), function(i)PCA$loadings[lookup==i,,drop=FALSE])
    names(loadings) <- names(X)
    mod <- list(blockLoadings=loadings, scores=scores, samplelinked=samplelinked)
    mod$info <- list(method = "Simultaneous Component Analysis", 
                     scores = "Common scores", loadings = "Not used",
                     blockScores = "Not used", blockLoadings = "Block loadings")
    attr(mod$scores, 'explvar') <- PCA$explvar
    for(i in 1:length(X)){
      explI <- numeric(length(PCA$explvar))
      denoI <- base::norm(X[[i]],'F')^2
      for(j in 1:length(PCA$explvar)){
        explI[j] <- sum(diag(scores[,1:j,drop=FALSE] %*% tcrossprod(crossprod(loadings[[i]][,1:j,drop=FALSE]), scores[,1:j,drop=FALSE])))/denoI*100
      }
      attr(mod$blockLoadings[[i]],'explvar') <- diff(c(0,explI))
    }
  }
  mod$explvar <- PCA$explvar
  mod$call <- match.call()
  mod$data <- list(X = X)
  class(mod) <- c('multiblock','list')
  return(mod)
}

#' Generalized Canonical Analysis - GCA
#' @param X \code{list} of input blocks.
#' @param ncomp \code{integer} number of components to extract, either single integer (equal for all blocks), vector (individual per block) or 'max' for maximum possible number of components.
#' @param svd \code{logical} indicating if Singular Value Decomposition approach should be used (default=TRUE).
#' @param tol \code{numeric} tolerance for component inclusion (singular values).
#' @param corrs \code{logical} indicating if correlations should be calculated for RGCCA based approach.
#' @param ... additional arguments for RGCCA approach.
#' 
#' @description This is an interface to both SVD-based (default) and RGCCA-based GCA (wrapping the
#' \code{RGCCA::rgcca} function)
#' 
#' @details GCA is a generalisation of Canonical Correlation Analysis to handle three or more
#' blocks. There are several ways to generalise, and two of these are available through \code{gca}.
#' The default is an SVD based approach estimating a common subspace and measuring mean squared
#' correlation to this. An alternative approach is available through RGCCA. For the SVD based
#' approach, the \code{ncomp} parameter controls the block-wise decomposition while the following
#' the consensus decomposition is limited to the minimum number of components from the individual blocks.
#' 
#' @return \code{multiblock} object including relevant scores and loadings. Relevant plotting functions: \code{\link{multiblock_plots}} 
#' and result functions: \code{\link{multiblock_results}}. \code{blockCoef} contains canonical coefficients, while
#' \code{blockDecomp} contains decompositions of each block.
#' 
#' @references
#' * Carroll, J. D. (1968). Generalization of canonical correlation analysis to three or more sets of variables. Proceedings of the American Psychological Association, pages 227-22.
#' * Van der Burg, E. and Dijksterhuis, G. (1996). Generalised canonical analysis of individual sensory profiles and instrument data, Elsevier, pp. 221–258.
#' 
#' @examples 
#' data(potato)
#' potList <- as.list(potato[c(1,2,9)])
#' pot.gca <- gca(potList)
#' plot(scores(pot.gca), labels="names")
#'
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results and plotting are found in \code{\link{multiblock_results}} and \code{\link{multiblock_plots}}, respectively.
#' @export
gca <- function(X, ncomp='max', svd=TRUE, tol=10^-12, corrs=TRUE, ...){
  if(length(ncomp)==1 && ncomp=='max'){
    ncomp <- unlist(lapply(X, function(x)min(nrow(x)-1, ncol(x))))
  }
  if(svd){
    obj <- gca.svd(X=X, tol=tol, ncomp=ncomp)
    obj$call = match.call()
    obj$data <- list(X = X)
    return(obj)
  } else {
    obj <- gca.rgcca(X=X, scale=FALSE, ncomp=ncomp, corrs=corrs, ...)
    obj$call = match.call()
    obj$data <- list(X = X)
    return(obj)
  }
}
gca.rgcca <- function(X, scale=FALSE, ncomp='max', corrs=TRUE, ...){
  n_block <- length(X)
  if(length(ncomp)==1){ ncomp <- rep(ncomp,n_block+1) }
  if(length(ncomp)==n_block){ ncomp <- c(ncomp, min(ncomp))}
  X <- lapply(X, function(i) scale(i, scale = FALSE))
  X[[n_block+1]] <- do.call(cbind, X)
  C <- matrix(0, n_block+1, n_block+1)
  C[n_block+1,1:n_block] <- 1; C[1:n_block,n_block+1] <- 1
  res <- RGCCA::rgcca(A = X, C = C, tau=rep(0,n_block+1), verbose = FALSE, scale = FALSE, ncomp=ncomp, scheme = "factorial", ...)
  # A = coefficients and global coefficients
  if(corrs){
    # Correlations: diag(cor(X[[i]]%*%res$astar[[i]],X[[j]]%*%res$astar[[j]]))
    lb <- unlist(lapply(1:(n_block-1), function(i)lapply((i+1):n_block, function(j)paste0(i,'-',j))))
    cs <- matrix(unlist(lapply(1:(n_block-1), function(i)lapply((i+1):n_block, function(j)diag(cor(X[[i]]%*%res$astar[[i]],X[[j]]%*%res$astar[[j]]))))), nrow=min(ncomp))
    colnames(cs) <- lb; rownames(cs) <- paste0(1:min(ncomp)," comp.")
    return(list(A = res$astar, B = NULL, corrs = cs, X = X, rgcca = res))
  } else {
    return(list(A = res$astar, B = NULL, X = X, rgcca = res))
  }
}
gca.svd <- function(X, tol=10^-12, ncomp=1){
  n_block <- length(X)
  if(length(ncomp)==1){ ncomp <- rep(ncomp,n_block) }
  n <- nrow(X[[1]])
  X <- lapply(X, function(i) scale(i, scale = FALSE))
  minrank <- min(min(unlist(lapply(X,ncol))),n)
  T <- Pb <- list()
  for(i in 1:n_block){
    udv <- svd(X[[i]])
    thisrank <- ifelse(sum(udv$d>tol) > 1, sum(udv$d>tol), 1)
    if(thisrank < ncomp[i]){
      warning(paste0("'ncomp' reduced due to low singular value for block ",i))
    } else {
      thisrank <- ncomp[i]
    }
    minrank <- min(thisrank, minrank)
    T[[i]]  <- udv$u[,1:thisrank,drop=FALSE]/rep(apply(udv$u[,1:thisrank,drop=FALSE],2,sd),each=n)
    Pb[[i]] <- udv$v[,1:thisrank,drop=FALSE] * rep(udv$d[1:thisrank], each=nrow(udv$v))
  }
  
  udv <- svd(do.call(cbind,T))
  C <- udv$u[,1:minrank,drop=FALSE]
  P <- udv$v[,1:minrank,drop=FALSE] * rep(udv$d[1:minrank], each=nrow(udv$v))
  A <- B <- U <- list()
  for(i in 1:n_block){
    A[[i]] <- pracma::pinv(X[[i]]) %*% C
    U[[i]] <- X[[i]] %*% A[[i]]  }
  
  ii <- 0
  R <- numeric(minrank)
  for(i in 1:(n_block-1)){
    for(j in (i+1):n_block){
      ii <- ii+1
      R <- R + diag(cor(U[[i]],U[[j]]))
    }
  }
  R <- R/ii
  info <- list(method = "Generalised Canonical Analysis", 
               scores = "Consensus scores", loadings = "Consensus loadings",
               blockScores = "Canonical scores", blockLoadings = "Not used")
  colnames(C) <- colnames(P) <- paste0('Comp ', 1:ncol(C))
  rownames(C) <- rownames(X[[1]])
  for(i in 1:length(X)){
    colnames(U[[i]])  <- paste0('Comp ', 1:ncol(U[[i]]))
    colnames(A[[i]]) <- paste0('Comp ', 1:ncol(A[[i]]))
    rownames(U[[i]])  <- rownames(X[[i]])
    rownames(A[[i]]) <- colnames(X[[i]])
  }
  names(U) <- names(A) <- names(X)
  rownames(P) <- paste0("Consensus dim ", 1:nrow(P))
  obj <- list(scores=C, loadings=P, blockScores=U, 
              blockCoef=A, blockDecomp=list(scores=T,loading=Pb), cor=R, info=info)
  class(obj) <- list('multiblock','list')
  return(obj)
}

#' Generalized Procrustes Analysis - GPA
#' @param X \code{list} of input blocks.
#' @param graph \code{logical} indicating if decomposition should be plotted.
#' @param ... additional arguments for RGCCA approach.
#' 
#' @description This is a wrapper for the \code{FactoMineR::GPA} function for computing GPA.
#' 
#' @details GPA is a generalisation of Procrustes analysis, where one matrix is scaled and 
#' rotated to be as similar as possible to another one. Through the generalisation, individual
#' scaling and rotation of each input matrix is performed against a common
#' representation which is estimated in an iterative manner.
#' 
#' @return \code{multiblock} object including relevant scores and loadings. Relevant plotting functions: \code{\link{multiblock_plots}} 
#' and result functions: \code{\link{multiblock_results}}.
#' 
#' @references Gower, J. C. (1975). Generalized procrustes analysis. Psychometrika. 40: 33–51.
#' 
#' @examples
#' data(potato)
#' potList <- as.list(potato[c(1,2,9)])
#' pot.gpa    <- gpa(potList)
#' plot(scores(pot.gpa), labels="names")
#'
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results and plotting are found in \code{\link{multiblock_results}} and \code{\link{multiblock_plots}}, respectively.
#' @export
gpa <- function(X, graph = FALSE, ...){
  Xcat <- as.data.frame(do.call(cbind,X))
  ret  <- FactoMineR::GPA(Xcat, group = unlist(lapply(X,ncol)), graph = graph, ...)
  blockScores <- blockLoadings <- list()
  for(i in 1:length(X)){
    usv <- svd(X[[i]] - rep(colMeans(X[[i]]), each=nrow(X[[i]])))
    blockScores[[i]]   <- usv$u * rep(usv$d, each=nrow(X[[i]]))
    blockLoadings[[i]] <- usv$v
    dimnames(blockScores[[i]])   <- list(rownames(X[[i]]), paste0('Comp ',1:ncol(blockScores[[i]])))
    dimnames(blockLoadings[[i]]) <- list(colnames(X[[i]]), paste0('Comp ',1:ncol(blockScores[[i]])))
    explvar <- 100*(usv$d^2/sum(usv$d^2)); names(explvar) <- paste0('Comp ',1:ncol(blockScores[[i]]))
    attr(blockScores[[i]], 'explvar') <- attr(blockLoadings[[i]], 'explvar') <- explvar
  }
  names(blockScores) <- names(blockLoadings) <- names(X)
  scores   <- ret$consensus
  loadings <- diag(ncol(ret$consensus))
  dimnames(scores)   <- list(rownames(X[[1]]), paste0('Comp ',1:ncol(scores)))
  dimnames(loadings) <- list(paste0('Dummy ',1:ncol(scores)), paste0('Comp ',1:ncol(scores)))
  info <- list(method = "Generalised Procrustes Analysis", 
               scores = "Consensus scores", loadings = "Not used",
               blockScores = "PCA scores", blockLoadings = "PCA loadings")
  obj <- list(scores = scores, loadings = loadings, 
              blockScores = blockScores, blockLoadings = blockLoadings, 
              info = info, GPA = ret, call = match.call(), data = list(X = X))
  class(obj) <- c('multiblock','list')
  return(obj)
}

#' Multiple Factor Analysis - MFA
#' @param X \code{list} of input blocks.
#' @param type \code{character} vector indicating block types, defaults to \code{rep("c", length(X))} for continuous values.
#' @param graph \code{logical} indicating if decomposition should be plotted.
#' @param ... additional arguments for RGCCA approach.
#' 
#' @description This is a wrapper for the \code{FactoMineR::MFA} function for computing MFA.
#' 
#' @details MFA is a methods typically used to compare several equally sized matrices. It is
#' often used in sensory analyses, where matrices consist of sensory characteristics and products,
#' and each assessor generates one matrix each. In its basic form, MFA scales all matrices by their
#' largest eigenvalue, concatenates them and performs PCA on the result. There are several
#' possibilities for plots and inspections of the model, handling of categorical and continuous
#' inputs etc. connected to MFA.
#' 
#' @return \code{multiblock} object including relevant scores and loadings. Relevant plotting functions: \code{\link{multiblock_plots}} 
#' and result functions: \code{\link{multiblock_results}}.
#' 
#' @references Pagès, J. (2005). Collection and analysis of perceived product inter-distances using multiple factor analysis: Application to the study of 10 white wines from the Loire valley. Food Quality and Preference, 16(7), 642–649.
#' 
#' @examples 
#' data(potato)
#' potList <- as.list(potato[c(1,2,9)])
#' pot.mfa    <- mfa(potList)
#' if(interactive()){
#'   plot(pot.mfa$MFA)
#' }
#'
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results and plotting are found in \code{\link{multiblock_results}} and \code{\link{multiblock_plots}}, respectively.
#' @export
mfa <- function(X, type = rep("c", length(X)), graph = FALSE, ...){
  ret <- FactoMineR::MFA(as.data.frame(do.call(cbind,X)), unlist(lapply(X,ncol)), type = type, graph = graph, ...)
  info <- list(method = "Multiple Factor Analysis", 
               scores = "Global scores", loadings = "Global loadings",
               blockScores = "Individual PCA scores", blockLoadings = "Individual PCA loadings")
  blockLoadings <- lapply(ret$separate.analyses, function(x)x$svd$V)
  blockScores   <- lapply(ret$separate.analyses, function(x)x$ind$coord)
  for(i in 1:length(X))
    rownames(blockLoadings[[i]]) <- colnames(X[[i]])
  names(blockLoadings) <- names(blockScores) <- names(X)
  obj <- list(scores = ret$ind$coord, loadings = ret$global.pca$svd$V,
              blockScores   = blockScores,
              blockLoadings = blockLoadings,
              info = info, MFA = ret, call = match.call())
  colnames(obj$scores) <- colnames(obj$loadings) <- paste0('Comp ', 1:ncol(obj$scores))
  rownames(obj$loadings) <- unlist(lapply(blockLoadings, rownames))
  obj$explvar <- attr(obj$scores, 'explvar') <- attr(obj$loadings, 'explvar') <- 100*(ret$global.pca$svd$vs^2/sum(ret$global.pca$svd$vs))
  for(i in 1:length(X)){
    explvar <- 100*(ret$separate.analyses[[i]]$svd$vs^2/sum(ret$separate.analyses[[i]]$svd$vs^2))
    attr(obj$blockScores[[i]], 'explvar') <- attr(obj$blockLoadings[[i]], 'explvar') <- explvar
    colnames(obj$blockScores[[i]]) <- colnames(obj$blockLoadings[[i]]) <- paste0('Comp ',1:ncol(obj$blockScores[[i]]))
  }
  obj$data <- list(X = X)
  class(obj) <- c('multiblock','list')
  return(obj)
}

#' PCA-GCA
#' @param X \code{list} of input blocks
#' @param commons \code{numeric} giving the highest number of blocks to combine when calculating local or common scores.
#' @param auto \code{logical} indicating if automatic choice of complexities should be used.
#' @param auto.par \code{named list} setting limits for automatic choice of complexities.
#' @param manual.par \code{named list} for manual choice of blocks. The list consists of \code{ncomp} which indicates the number of components to extract from each block and \code{ncommon} which is the corresponding for choosing the block combinations (local/common). For the latter, use unique_combos(n_blocks, commons) to see order of local/common blocks. Component numbers will be reduced if simpler models give better predictions. See example.
#' @param tol \code{numeric} tolerance for component inclusion (singular values).
#' 
#' @return \code{multiblock} object including relevant scores and loadings. Relevant plotting functions: \code{\link{multiblock_plots}} 
#' and result functions: \code{\link{multiblock_results}}. Distinct components are marked as 'D(x), Comp c' for block x and component c
#' while local and common components are marked as "C(x1, x2), Comp c", where x1 and x2 (and more) are block numbers. 
#' 
#' @description PCA-GCA is a methods which aims at estimating subspaces of common, local and
#' distinct variation from two or more blocks. 
#' 
#' @details The name PCA-GCA comes from the process of first applying
#' PCA to each block, then using GCA to estimate local and common components, and finally
#' orthogonalising the block-wise scores on the local/common ones and re-estimating these to
#' obtain distinct components. The procedure is highly similar to the supervised method
#' PO-PLS, where the PCA steps are exchanged with PLS.
#' 
#' @references Smilde, A., Måge, I., Naes, T., Hankemeier, T.,Lips, M., Kiers, H., Acar, E., and Bro, R.(2017). Common and distinct components in data fusion. Journal of Chemometrics, 31(7), e2900.
#' 
#' @examples 
#' data(potato)
#' potList <- as.list(potato[c(1,2,9)])
#' pot.pcagca <- pcagca(potList)
#' 
#' # Show origin and type of all components
#' lapply(pot.pcagca$blockScores,colnames)
#' 
#' # Basic multiblock plot
#' plot(scores(pot.pcagca, block=2), comps=1, labels="names")
#'
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results and plotting are found in \code{\link{multiblock_results}} and \code{\link{multiblock_plots}}, respectively.
#' @export
pcagca <- function(X, commons=2, auto=TRUE, auto.par=list(explVarLim=40, rLim=0.8),
                   manual.par=list(ncomp=0, ncommon=0), tol=10^-12){
  
  # commons can be a list of integer vectors specifying which local/common block combinations should be included or a single integer indicating the highest order of commonness between blocks, e.g., 2 means all combinations of two blocks are included.
  
  # Initialize
  n_block <- length(X)
  n <- nrow(X[[1]])
  X <- lapply(X, scale, scale=FALSE)
  totVar <- lapply(X, function(x)sum(x^2))
  Scores <- Loadings <- explVar <- U <- T <- UC <- blockLabels <- list()
  
  # Check components
  if(auto){
    ncomp <- pmin(unlist(lapply(X,ncol)),nrow(X[[1]])-1)
    ncommon <- list()
  } else {
    ncomp <- manual.par$ncomp
    if(length(ncomp)==1){
      ncomp <- rep(ncomp, n_block)
    }
    ncommon <- manual.par$ncommon
  }
  
  # List of common block combinations
  if(is.numeric(commons)){
    commons <- unique_combos(n_block, commons)
  }
  commonLabels <- lapply(commons, function(x)paste(x,sep="",collapse=","))
  
  # Prepare for storage
  for(i in 1:n_block){
    Scores[[i]] <- matrix(0, n,0)
    Loadings[[i]] <- matrix(0, ncol(X[[i]]),0)
  }
  
  # Basis scores for each block
  for(i in 1:n_block){
    udv <- svd(X[[i]], ncomp[i], ncomp[i])
    # U[[i]] <- udv$u[,1:ncomp[i],drop=FALSE]
    T[[i]] <- udv$u[,1:ncomp[i],drop=FALSE] %*% diag(udv$d[1:ncomp[i]])
    ncomp[i] <- pcaopt(X[[i]], T[[i]], udv$v, ncomp[i])
    T[[i]] <- T[[i]][,1:ncomp[i],drop=FALSE]
    blockLabels[[i]] <- character(0)
  }
  
  # Common components from GCA
  n_common <- length(commons)
  R <- numeric(0)
  for(i in 1:n_common){
    gcaComp <- gca(T[commons[[i]]])
    UC[[i]] <- gcaComp$blockScores
    
    # Explained variance
    explVarX <- matrix(0, length(gcaComp$cor),length(commons[[i]]))
    for(j in 1:length(commons[[i]])){
      XX <- X[[commons[[i]][j]]]
      for(k in 1:length(gcaComp$cor)){
        ss <- gcaComp$blockScores[[j]][,k]/c(sqrt(crossprod(gcaComp$blockScores[[j]][,k])))
        loads <- crossprod(XX, ss)
        XX <- XX - tcrossprod(ss, loads)
        explVarX[k,j] <- (totVar[[commons[[i]][j]]]-sum(XX^2))/totVar[[commons[[i]][j]]]*100
      }
    }
    explVarX <- diff(rbind(numeric(length(commons[[i]])), explVarX))
    if(auto){
      nc <- which(gcaComp$cor>auto.par$rLim & rowMeans(explVarX)>auto.par$explVarLim )
      ncommon[[i]] <- nc
    }

    if((!length(ncommon[[i]])==0) && ncommon[[i]]>0){
      R <- c(R, gcaComp$cor[ncommon[[i]]])
      for(j in 1:length(commons[[i]])){
        # Store common scores
        u <- gcaComp$blockScores[[j]][,ncommon[[i]],drop=FALSE]
        Scores[[commons[[i]][j]]] <- cbind(Scores[[commons[[i]][j]]], u/rep(sqrt(diag(crossprod(u))), each=n))
        for(k in 1:length(ncommon[[i]])){
          blockLabels[[commons[[i]][j]]] <- c(blockLabels[[commons[[i]][j]]], paste("C(", commonLabels[[i]],"), Comp ",k,sep=""))
        }

        # Deflate basis scores
        T[[commons[[i]][j]]] <- T[[commons[[i]][j]]] - u%*%solve(crossprod(u))%*%crossprod(u,T[[commons[[i]][j]]])
      }
    }
  }
  
  # Recompose basis into unique scores
  for(i in 1:n_block){
    udv <- svd(T[[i]])
    Scores[[i]] <- cbind(Scores[[i]], udv$u[,udv$d>tol, drop=FALSE])
    for(k in 1:sum(udv$d>tol)){
      blockLabels[[i]] <- c(blockLabels[[i]], paste("D(",i,"), Comp ",k,sep=""))
    }
  }
  
  # Calculate loadings
  for(i in 1:n_block){
    for(j in 1:ncol(Scores[[i]])){
      Loadings[[i]] <- cbind(Loadings[[i]],crossprod(X[[i]], Scores[[i]][,j,drop=FALSE]))
      X[[i]] <- X[[i]] - tcrossprod(Scores[[i]][,j,drop=FALSE], Loadings[[i]][,j,drop=FALSE])
      if(length(explVar)<i){
        explVar[[i]] <- numeric(0)}
      explVar[[i]] <- c(explVar[[i]], (totVar[[i]]-sum(X[[i]]^2))/totVar[[i]]*100)
    }
  }
  names(ncommon) <- commonLabels
  for(i in 1:n_block){
    rownames(Scores[[i]]) <- rownames(X[[i]])
  }
  explVar  <- lapply(1:n_block, function(i){e <- diff(c(0,explVar[[i]])); names(e) <- blockLabels[[i]];e})
  Scores   <- lapply(1:n_block, function(i){e <- Scores[[i]]; colnames(e) <- blockLabels[[i]];e})
  Loadings <- lapply(1:n_block, function(i){e <- Loadings[[i]]; colnames(e) <- blockLabels[[i]];e})
  if(auto){
    Scores   <- lapply(1:n_block,function(i)Scores[[i]][,explVar[[i]]>auto.par$explVarLim, drop=FALSE])
    Loadings <- lapply(1:n_block,function(i)Loadings[[i]][,explVar[[i]]>auto.par$explVarLim, drop=FALSE])
    explVar  <- lapply(1:n_block,function(i)explVar[[i]][explVar[[i]]>auto.par$explVarLim, drop=FALSE])
  }
  names(Scores) <- names(Loadings) <- names(explVar) <- names(X)
  for(i in 1:n_block){
    attr(Scores[[i]], 'explvar') <- attr(Loadings[[i]], 'explvar') <- explVar[[i]]
  }
  info <- list(method = "PCA-GCA", 
               scores = "Not available", loadings = "Not available",
               blockScores = "Common and distinct scores", blockLoadings = "Common and distinct loadings")
  obj <- list(blockScores=Scores, blockLoadings=Loadings, R=R, explVar=explVar, 
              ncomp=unlist(lapply(explVar, length)), ncommon=ncommon, 
              info=info, call=match.call(), data = list(X = X))
  class(obj) <- c("multiblock","list")
  return(obj)
}

#' Distinctive and Common Components with SCA - DISCO
#' @param X \code{list} of input blocks.
#' @param ncomp \code{integer} number of components to extract.
#' @param ... additional arguments (not used).
#' 
#' @return \code{multiblock} object including relevant scores and loadings. Relevant plotting functions: \code{\link{multiblock_plots}} 
#' and result functions: \code{\link{multiblock_results}}.
#' 
#' @description This is a wrapper for the \code{RegularizedSCA::DISCOsca} function for computing DISCO.
#' 
#' @details DISCO is a restriction of SCA where Alternating Least Squares is used for
#' estimation of loadings and scores. The SCA solution is rotated towards loadings (in sample linked mode) which are filled with 
#' zeros in a pattern resembling distinct, local and common components.
#' When used in sample linked mode and only selecting distinct components, it shares a 
#' resemblance to SO-PLS, only in an unsupervised setting. Explained variances
#' are computed as proportion of block varation explained by scores*loadings'.
#' 
#' @references Schouteden, M., Van Deun, K., Wilderjans, T. F., & Van Mechelen, I. (2014). Performing DISCO-SCA to search for distinctive and common information in linked data. Behavior research methods, 46(2), 576-587.
#' 
#' @examples
#' data(potato)
#' potList <- as.list(potato[c(1,2,9)])
#' pot.disco  <- disco(potList)
#' plot(scores(pot.disco), labels="names")
#'
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' @export
disco <- function(X, ncomp = 2, ...){
  Xc <- do.call(cbind,X)
  ret <- RegularizedSCA::DISCOsca(Xc, ncomp, unlist(lapply(X,ncol)))
  compNames <- character(ncomp)
  for(i in 1:ncomp){
    if(sum(ret$comdist[[1]][,i]) == 1)
      compNames[i] <- paste0('Comp ', i, ', D(', which(ret$comdist[[1]][,i]==1), ')')
    else
      compNames[i] <- paste0('Comp ', i, ', C(', paste(which(ret$comdist[[1]][,i]==1), collapse=",", sep=""), ')')
  }
  blockLoadings <- list(); j <- 0
  for(i in 1:length(X)){
    blockLoadings[[i]] <- ret$Prot_best[[1]][j+(1:ncol(X[[i]])),,drop=FALSE]
    j <- j+ncol(X[[i]])
  }
  colnames(ret$Trot_best[[1]]) <- colnames(ret$Prot_best[[1]]) <- compNames
  for(i in 1:length(X)){
    colnames(blockLoadings[[i]]) <- compNames
    rownames(blockLoadings[[i]]) <- colnames(X[[i]])
  }
  rownames(ret$Trot_best[[1]]) <- rownames(X[[1]])
  rownames(ret$Prot_best[[1]]) <- unlist(lapply(blockLoadings, rownames))
  names(blockLoadings) <- names(X)
  info <- list(method = "Distinctive and Common Components with SCA",
               scores = "Scores", loadings = "Concatenated loadings",
               blockScores = "Not used", blockLoadings = "Block-wise loadings")
  obj <- list(scores = ret$Trot_best[[1]], loadings = ret$Prot_best[[1]], blockLoadings = blockLoadings,
              info = info, DISCOsca = ret, call = match.call(), data = list(X = X))
  xFro <- unlist(lapply(X, function(x)base::norm(scale(x,scale=FALSE),type='F')^2))
  explvar <- obj$DISCOsca$propExp_component[[1]]
  for(j in 1:dim(explvar)[1]){
    x <- scale(X[[j]], scale=FALSE)
    for(i in 1:ncomp){
      explvar[j,i] <- 100-norm(x-obj$scores[,i,drop=FALSE]%*%t(obj$blockLoadings[[j]][,i,drop=FALSE]),type="F")^2/xFro[j]*100
    }
  }
  dimnames(explvar) <- list(names(X),colnames(obj$loadings))
  for(i in 1:length(X)){
    attr(obj$blockLoadings[[i]], 'explvar') <- explvar[i,] # diff(c(0,obj$DISCOsca$propExp_component[[1]][i,]))/xFro[[i]]*100
  }
  obj$explvar <- attr(obj$scores, "explvar") <- attr(obj$loadings, "explvar") <- explvar#diff(c(0,diag(crossprod(obj$loadings))))/sum(xFro)*100
  class(obj) <- c("multiblock","list")
  return(obj)
}

#' Hierarchical Principal component analysis - HPCA
#' @param X \code{list} of input blocks.
#' @param ncomp \code{integer} number of components to extract.
#' @param scale \code{logical} indicating if variables should be scaled.
#' @param verbose \code{logical} indicating if diagnostic information should be printed.
#' @param ... additional arguments for RGCCA.
#' 
#' @description This is a wrapper for the \code{RGCCA::rgcca} function for computing HPCA.
#' 
#' @details HPCA is a hierarchical PCA analysis which combines two or more blocks
#' into a two-level decomposition with block-wise loadings and scores and superlevel
#' common loadings and scores. The method is closely related to the supervised method MB-PLS
#' in structure.
#' 
#' @return \code{multiblock} object including relevant scores and loadings. Relevant plotting functions: \code{\link{multiblock_plots}} 
#' and result functions: \code{\link{multiblock_results}}.
#'
#' @references Westerhuis, J.A., Kourti, T., and MacGregor,J.F. (1998). Analysis of multiblock and hierarchical PCA and PLS models. Journal of Chemometrics, 12, 301–321.
#' 
#' @examples
#' data(potato)
#' potList <- as.list(potato[c(1,2,9)])
#' pot.hpca   <- hpca(potList)
#' plot(scores(pot.hpca), labels="names")
#'
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results and plotting are found in \code{\link{multiblock_results}} and \code{\link{multiblock_plots}}, respectively.
#' @export
hpca <- function(X, ncomp=2, scale=FALSE, verbose=FALSE, ...){
  n_block <- length(X)
  if(length(ncomp)==1){ ncomp <- rep(ncomp,n_block+1) }
  X <- lapply(X, function(i) scale(i, scale = FALSE))
  X[[n_block+1]] <- do.call(cbind, X)
  C <- matrix(0, n_block+1, n_block+1)
  C[n_block+1,1:n_block] <- 1; C[1:n_block,n_block+1] <- 1
  res <- RGCCA::rgcca(A = X, C = C, tau=c(rep(1,n_block),0), scale = scale, ncomp=ncomp, scheme = function(x) x^4, verbose = verbose, ...)
  # A = scores and superscore (last), B = loadings and superloadings (last)
  info <- list(method = "Hierarchical PCA",
               scores = "Super scores", loadings = "Super loadings",
               blockScores = "Block scores", blockLoadings = "Block loadings")
  scores   <- X[[n_block+1]]%*%res$astar[[n_block+1]]
  loadings <- res$astar[[n_block+1]]
  colnames(scores) <- colnames(loadings) <- paste0('Comp ', 1:ncomp[1])
  blockScores   <- colnamesList(lapply(1:n_block, function(i)X[[i]]%*%res$astar[[i]]), paste0('Comp ', 1:ncomp[1]))
  blockLoadings <- colnamesList(res$astar[1:n_block], paste0('Comp ', 1:ncomp[1]))
  for(i in 1:n_block)
    rownames(blockLoadings[[i]]) <- colnames(X[[i]])
  names(blockScores) <- names(blockLoadings) <- names(X[1:n_block])
  obj <- list(scores = scores, loadings = loadings, 
              blockScores = blockScores, blockLoadings = blockLoadings, X = X, 
              info = info, rgcca = res, call = match.call())
  for(i in 1:n_block){
    attr(obj$blockScores[[i]], "explvar") <- attr(obj$blockLoadings[[i]], "explvar") <- obj$rgcca$AVE$AVE_X[[i]]
  }
  obj$explvar <- res$AVE_outer_model
  obj$data <- list(X = X)
  class(obj) <- c("multiblock","list")
  return(obj)
}

#' Multiple Co-Inertia Analysis - MCOA
#' @param X \code{list} of input blocks.
#' @param ncomp \code{integer} number of components to extract.
#' @param scale \code{logical} indicating if variables should be scaled.
#' @param verbose \code{logical} indicating if diagnostic information should be printed.
#' @param ... additional arguments for RGCCA.
#' 
#' @return \code{multiblock} object including relevant scores and loadings. Relevant plotting functions: \code{\link{multiblock_plots}} 
#' and result functions: \code{\link{multiblock_results}}.
#' 
#' @description This is a wrapper for the \code{RGCCA::rgcca} function for computing MCOA.
#' 
#' @details MCOA resembles GCA and MFA in that it creates a set of reference scores, for which each
#' block's individual scores should correlate maximally too, but also the variance within
#' each block should be taken into account. A single component solution is equivalent to a
#' PCA on concatenated blocks scaled by the so called inverse inertia.
#' 
#' @references 
#' * Le Roux; B. and H. Rouanet (2004). Geometric Data Analysis, From Correspondence Analysis to Structured Data Analysis. Dordrecht. Kluwer: p.180.
#' * Greenacre, Michael and Blasius, Jörg (editors) (2006). Multiple Correspondence Analysis and Related Methods. London: Chapman & Hall/CRC.
#' 
#' @examples 
#' data(potato)
#' potList <- as.list(potato[c(1,2,9)])
#' pot.mcoa   <- mcoa(potList)
#' plot(scores(pot.mcoa), labels="names")
#'
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results and plotting are found in \code{\link{multiblock_results}} and \code{\link{multiblock_plots}}, respectively.
#' @export
mcoa <- function(X, ncomp=2, scale=FALSE, verbose=FALSE, ...){
  n_block <- length(X)
  if(length(ncomp)==1){ ncomp <- rep(ncomp,n_block+1) }
  X <- lapply(X, function(i) scale(i, scale = FALSE))
  X[[n_block+1]] <- do.call(cbind, X)
  C <- matrix(0, n_block+1, n_block+1)
  C[n_block+1,1:n_block] <- 1; C[1:n_block,n_block+1] <- 1
  res <- RGCCA::rgcca(A = X, C = C, tau=c(rep(1,n_block),0), verbose = verbose, scale = scale, ncomp=ncomp, scheme = "factorial", ...)
  # A = coefficients and global coefficients
  # return(list(A = res$astar, B = NULL, X = X, rgcca = res))
  info <- list(method = "Multiple Co-Inertia Analysis",
               scores = "Super scores", loadings = "Super loadings",
               blockScores = "Block scores", blockLoadings = "Block loadings")
  scores   <- X[[n_block+1]]%*%res$astar[[n_block+1]]
  loadings <- res$astar[[n_block+1]]
  colnames(scores) <- colnames(loadings) <- paste0('Comp ', 1:ncomp[1])
  blockScores   <- colnamesList(lapply(1:n_block, function(i)X[[i]]%*%res$astar[[i]]), paste0('Comp ', 1:ncomp[1]))
  blockLoadings <- colnamesList(res$astar[1:n_block], paste0('Comp ', 1:ncomp[1]))
  for(i in 1:n_block)
    rownames(blockLoadings[[i]]) <- colnames(X[[i]])
  names(blockScores) <- names(blockLoadings) <- names(X[1:n_block])
  obj <- list(scores = scores, loadings = loadings, 
              blockScores = blockScores, blockLoadings = blockLoadings, X = X, 
              info = info, rgcca = res, call = match.call())
  for(i in 1:n_block){
    attr(obj$blockScores[[i]], "explvar") <- attr(obj$blockLoadings[[i]], "explvar") <- obj$rgcca$AVE$AVE_X[[i]]
  }
  obj$explvar <- res$AVE_outer_model
  obj$data <- list(X = X)
  class(obj) <- c("multiblock","list")
  return(obj)
}

#' Joint and Individual Variation Explained - JIVE
#' @param X \code{list} of input blocks.
#' @param ... additional arguments for \code{r.jive::jive}.
#' 
#' @return \code{multiblock} object including relevant scores and loadings. Relevant plotting functions: \code{\link{multiblock_plots}} 
#' and result functions: \code{\link{multiblock_results}}.
#' 
#' @description This is a wrapper for the \code{r.jive::jive} function for computing JIVE.
#' 
#' @details Jive performs a decomposition of the variation in two or more blocks into
#' low-dimensional representations of individual and joint variation plus residual variation.
#' 
#' @references Lock, E., Hoadley, K., Marron, J., and Nobel, A. (2013) Joint and individual variation explained (JIVE) for integrated analysis of multiple data types. Ann Appl Stat, 7 (1), 523–542.
#' 
#' @examples 
#' \donttest{ # Too time consuming for testing
#' data(candies)
#' candyList <- lapply(1:nlevels(candies$candy),function(x)candies$assessment[candies$candy==x,])
#' can.jive  <- jive(candyList)
#' summary(can.jive)
#' }
#' 
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' @export
jive <- function(X, ...){
  r.jive::jive(X, ...)
}

#' Structuration des Tableaux à Trois Indices de la Statistique - STATIS
#' @param X \code{list} of input blocks.
#' @param ncomp \code{integer} number of components to extract.
#' @param scannf \code{logical} indicating if eigenvalue bar plot shoulde be displayed.
#' @param tol \code{numeric} eigenvalue threshold tolerance.
#' @param ... additional arguments (not used).
#' 
#' @return \code{multiblock} object including relevant scores and loadings. Relevant plotting functions: \code{\link{multiblock_plots}} 
#' and result functions: \code{\link{multiblock_results}}.
#' 
#' @description This is a wrapper for the \code{ade4::statis} function for computing STATIS.
#' 
#' @details STATIS is a method, related to MFA, for analysing two or more blocks. It also
#' decomposes the data into a low-dimensional subspace but uses a different scaling of the
#' individual blocks.
#' 
#' @references Lavit, C.; Escoufier, Y.; Sabatier, R.; Traissac, P. (1994). The ACT (STATIS method). Computational Statistics & Data Analysis. 18: 97
#' 
#' @examples
#' data(candies)
#' candyList <- lapply(1:nlevels(candies$candy),function(x)candies$assessment[candies$candy==x,])
#' can.statis <- statis(candyList)
#' plot(scores(can.statis), labels="names")
#'
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results and plotting are found in \code{\link{multiblock_results}} and \code{\link{multiblock_plots}}, respectively.
#' @export
statis <- function(X, ncomp = 3, scannf = FALSE, tol = 1e-07, ...){
  X_frame  <- as.data.frame(do.call(rbind, X))
  X_factor <- factor(unlist(lapply(1:length(X), function(x)rep(x,nrow(X[[x]])))))
  kta <- ktab.within(withinpca(X_frame, X_factor, scannf=scannf, nf=ncomp))
  ret <- ade4::statis(kta, scannf=scannf, nf=ncomp, tol=tol)
  scores <- as.matrix(ret$C.Co); loadings <- as.matrix(ret$C.li)
  blockScores <- list(); j <- 0
  for(i in 1:length(X)){
    blockScores[[i]] <- scores[j+(1:nrow(X[[i]])),,drop=FALSE]; j <- j+nrow(X[[i]])
  }
  names(blockScores) <- names(X)
  colnames(scores) <- colnames(loadings) <- paste0('Comp ', 1:ncomp)
  blockScores <- colnamesList(blockScores, paste0('Comp ', 1:ncomp))
  info <- list(method = "STATIS",
               scores = "Concatenated scores", loadings = "Loadings",
               blockScores = "Block-wise scores", blockLoadings = "Not used")
  obj <- list(scores = scores, loadings = loadings, blockScores = blockScores,
              info = info, statis = ret, call = match.call())
  obj$data <- list(X = X)
  class(obj) <- c("multiblock","list")
  return(obj)
  # Has plot and print in ade4
  # Clustatis in https://cran.r-project.org/web/packages/ClustBlock/ClustBlock.pdf
}

#' Higher Order Generalized SVD - HOGSVD
#' @param X \code{list} of input blocks.
#' 
#' @return \code{multiblock} object including relevant scores and loadings. Relevant plotting functions: \code{\link{multiblock_plots}} 
#' and result functions: \code{\link{multiblock_results}}.
#' 
#' @description This is a simple implementation for computing HOGSVD
#' 
#' @details HOGSVD is a generalisation of SVD to two or more blocks. It finds a common set
#' of loadings across blocks and individual sets of scores per block.
#' 
#' @references Ponnapalli, S. P., Saunders, M. A., Van Loan, C. F., & Alter, O. (2011). A higher-order generalized singular value decomposition for comparison of global mRNA expression from multiple organisms. PloS one, 6(12), e28072.
#' 
#' @examples
#' data(candies)
#' candyList <- lapply(1:nlevels(candies$candy),function(x)candies$assessment[candies$candy==x,])
#' can.hogsvd <- hogsvd(candyList)
#' scoreplot(can.hogsvd, block=1, labels="names")
#' 
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results and plotting are found in \code{\link{multiblock_results}} and \code{\link{multiblock_plots}}, respectively.
#' @export
hogsvd  <- function(X){
  # Assumes equal number of variables
  nvar  <- ncol(X[[1]])
  nblock <- length(X)
  
  # Calculate S
  A <- lapply(X, crossprod)
  Ainv <- lapply(A, solve)
  S <- matrix(0, nvar, nvar)
  for(i in 1:(nblock-1)){
    for(j in (i+1):nblock){
      S <- S + crossprod(A[[i]], Ainv[[j]]) + crossprod(A[[j]], Ainv[[i]])
    }
  }
  S <- S / (nblock * (nblock-1))
  
  # Eigen-decompose S
  eig <- eigen(S)
  eigen_values <- eig$values
  V <- eig$vectors
  
  # Calculate B
  Vinv <- solve(V)
  B <- lapply(X, function(x) t(tcrossprod(Vinv,x)))
  
  # Calculate U and sigmas
  vecnorm <- function(x)sqrt(crossprod(x))
  sigmas <- lapply(B, function(b)apply(b,2,vecnorm))
  U <- lapply(1:length(X), function(i)B[[i]]/rep(sigmas[[i]], each=nrow(X[[i]])))
  blockScores <-B
  names(blockScores) <- names(X)
  dimnames(V) <- list(colnames(X[[1]]), paste0('Comp ', 1:ncol(V)))
  for(i in 1:length(X))
    rownames(blockScores[[i]]) <- rownames(X[[i]])
  blockScores <- colnamesList(blockScores, paste0('Comp ', 1:ncol(V)))
  names(blockScores) <- names(X)
  info <- list(method = "Higher Order Generalised SVD",
               scores = "Not used", loadings = "Loadings",
               blockScores = "Block scores", blockLoadings = "Not used")
  obj <- list(loadings=V, blockScores=blockScores, bSnorm1=U, sigmas=sigmas, eigen_values=eigen_values,
              info = info, call = match.call())
  obj$explvar <- eigen_values^2/sum(eigen_values^2)*100
  obj$data <- list(X = X)
  class(obj) <- c("multiblock","list")
  return(obj)
}
