#' @title Parallel and Orthogonalised Partial Least Squares - PO-PLS
#' 
#' @description This is a basic implementation of PO-PLS with manual and automatic component selections.
#' 
#' @details PO-PLS decomposes a set of input data blocks into common, local and distinct components
#' through a process involving \code{\link{pls}} and \code{\link{gca}}. The \code{rLim} parameter is 
#' a lower bound for the GCA correlation when building common components, while explVarLim is the minimum
#' explained variance for common components and unique components.
#' 
#' @param X \code{list} of input blocks
#' @param Y \code{matrix} of response variable(s)
#' @param commons \code{numeric} giving the highest number of blocks to combine when calculating local or common scores.
#' @param auto \code{logical} indicating if automatic choice of complexities should be used.
#' @param auto.par \code{named list} setting limits for automatic choice of complexities. See Details.
#' @param manual.par \code{named list} for manual choice of blocks. The list consists of \code{ncomp} which indicates the number of components to extract from each block and \code{ncommon} which is the corresponding for choosing the block combinations (local/common). For the latter, use unique_combos(n_blocks, commons) to see order of local/common blocks. Component numbers will be reduced if simpler models give better predictions. See example.
#' 
#' @return A \code{multiblock} object with block-wise, local and common loadings and scores. Relevant plotting functions: \code{\link{multiblock_plots}} 
#' and result functions: \code{\link{multiblock_results}}.
#' 
#' @references 
#' * I Måge, BH Mevik, T Næs. (2008). Regression models with process variables and parallel blocks of raw material measurements. Journal of Chemometrics: A Journal of the Chemometrics Society 22 (8), 443-456
#' * I Måge, E Menichelli, T Næs (2012). Preference mapping by PO-PLS: Separating common and unique information in several data blocks. Food quality and preference 24 (1), 8-16
#'
#' @examples 
#' data(potato)
#' 
#' # Automatic analysis
#' pot.po.auto <- popls(potato[1:3], potato[['Sensory']][,1])
#' pot.po.auto$explVar
#' 
#' # Manual choice of up to 5 components for each block and 1, 0, and 2 blocks,
#' # respectively from the (1,2), (1,3) and (2,3) combinations of blocks.
#' pot.po.man <- popls(potato[1:3], potato[['Sensory']][,1], auto=FALSE, 
#'                 manual.par = list(ncomp=c(5,5,5), ncommon=c(1,0,2)))
#' pot.po.man$explVar
#' 
#' # Score plot for local (2,3) components
#' plot(scores(pot.po.man,3), comps=1:2, labels="names")
#' 
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results and plotting are found in \code{\link{multiblock_results}} and \code{\link{multiblock_plots}}, respectively.
#' @export
popls <- function(X, Y, commons=2, auto=TRUE, auto.par=list(explVarLim=40, rLim=0.8),
                   manual.par=list(ncomp=rep(0,length(X)), ncommon=list())){
  
  # commons can be a list of integer vectors specifying which local/common block combinations should be included or a single integer indicating the highest order of commonness between blocks, e.g., 2 means all combinations of two blocks are included.
  
  # Initialize
  n_block <- length(X)
  n <- nrow(X[[1]])
  X <- lapply(X, scale, scale=FALSE)
  totVar <- lapply(X, function(x)sum(x^2))
  Scores <- Loadings <- explVar <- U <- T <- UC <- blockLabels <- list()
  Y.res <- Y <- scale(as.matrix(Y), scale=FALSE)
  
  # Check components
  if(auto){
    ncomp <- pmin(unlist(lapply(X,ncol)),nrow(X[[1]])-2)
    ncommon <- list()
  } else {
    ncomp <- manual.par$ncomp
    if(length(ncomp)==1){
      ncomp <- rep(ncomp, n_block)
    }
    ncommon <- list()
    for(i in 1:length(manual.par$ncommon)){
      if(manual.par$ncommon[i]>0)
        ncommon[[i]] <- 1:manual.par$ncommon[i]
      else
        ncommon[[i]] <- 0
    }
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
    x <- X[[i]]
    pl <- plsr(Y ~ x, ncomp = max(ncomp[i],1), validation = "LOO")
    T[[i]] <- scores(pl)
    ncomp[i] <- which.min(apply(RMSEP(pl)$val,3,mean))-1
    # udv <- svd(X[[i]], ncomp[i], ncomp[i])
    # T[[i]] <- udv$u[,1:ncomp[i],drop=FALSE] %*% diag(udv$d[1:ncomp[i]])
    # ncomp[i] <- pcaopt(X[[i]], T[[i]], udv$v, ncomp[i])
    if(ncomp[i]>0){
      T[[i]] <- T[[i]][,1:ncomp[i],drop=FALSE]
    } else {
      T[[i]] <- c(0)
    }
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
    
    if(length(ncommon[[i]])>0 && any(ncommon[[i]]>0)){
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
        
        # Deflate response
        Y.res <- Y.res - u%*%solve(crossprod(u))%*%crossprod(u,Y.res)
      }
    }
  }
  
  # Recompose basis into unique scores
  for(i in 1:n_block){
    x <- T[[i]]
    pl <- plsr(Y.res ~ x, ncomp = ncomp[i], validation = "LOO")
    nc <- which.min(apply(RMSEP(pl)$val,3,mean))-1
    if(nc>0){
      uu <- scores(pl)[,1:nc,drop=FALSE]
      Scores[[i]] <- cbind(Scores[[i]], uu/rep(sqrt(diag(crossprod(uu))), each=n))
    }

    if(nc > 0){
      for(k in 1:nc){
        blockLabels[[i]] <- c(blockLabels[[i]], paste("D(",i,"), Comp ",k,sep=""))
      }
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
  info <- list(method = "PO-PLS", 
               scores = "Not available", loadings = "Not available",
               blockScores = "Common and distinct scores", blockLoadings = "Common and distinct loadings")
  obj <- list(blockScores=Scores, blockLoadings=Loadings, R=R, explVar=explVar, 
              ncomp=unlist(lapply(explVar, length)), ncommon=ncommon, 
              info=info, call=match.call(), opt.comp=ncomp)
  class(obj) <- c("multiblock","list")
  return(obj)
}
