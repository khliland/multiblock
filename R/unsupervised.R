#' @name unsupervised
#' @title Unsupervised Multiblock Methods
#' @aliases sca gca mfa pcagca disco jive statis hogsvd hpca mcoa
#' @importFrom RGCCA rgcca
#' @importFrom FactoMineR MFA GPA
#' @importFrom RegularizedSCA DISCOsca
#' @importFrom ade4 statis ktab.within withinpca
#' @importFrom r.jive jive summary.jive plot.jive
#' 
#' @description Collection of unsupervised multiblock methods:
#' * SCA - Simultaneous Component Analysis (\code{sca})
#' * GCA - Generalized Canonical Analysis (\code{gca})
#' * GPA - Generalized Procrustes Analysis (\code{gpa})
#' * MFA - Multiple Factor Analysis (\code{mfa})
#' * PCA-GCA (\code{pcagca})
#' * DISCO - Distinctive and Common Components with SCA (\code{disco})
#' * HPCA - Hierarchical Principal component analysis (\code{hpca})
#' * MCOA - Multiple Co-Inertia Analysis (\code{mcoa})
#' * JIVE - Joint and Individual Variation Explained (\code{jive})
#' * STATIS - Structuration des Tableaux à Trois Indices de la Statistique (\code{statis})
#' * HOGSVD - Higher Order Generalized SVD (\code{hogsvd})
#' 
#' @details 
#' Original documentation of STATIS: \link[ade4]{statis}.
#' JIVE, STATIS and HOGSVD assume variable linked matrices/data.frames, while SCA handles both links.
#' 
#' @examples 
#' # Object linked data
#' potList <- as.list(potato[c(1,2,9)])
#' pot.sca    <- sca(potList)
#' pot.gca    <- gca(potList)
#' pot.gpa    <- gpa(potList)
#' pot.mfa    <- mfa(potList)
#' pot.pcagca <- pcagca(potList)
#' pot.disco  <- disco(potList)
#' pot.hpca   <- hpca(potList)
#' pot.mcoa   <- mcoa(potList)
#' 
#' # Variable linked data
#' candyList <- lapply(1:nlevels(candies$candy),function(x)candies$assessment[candies$candy==x,])
#' can.jive   <- jive(candyList)
#' can.statis <- statis(candyList)
#' can.hogsvd <- hogsvd(candyList)
#' plot(can.statis$statis)
#
# Include also Penalized Exponential SCA (SLIDE, penalty-based)?????????
# 
#' @rdname unsupervised
#' @export
sca <- function(X, ncomp=2, scale=FALSE, samplelinked = 'auto', ...){
  # SVD/PCA based SVD-P with block-centering (and global scaling)
  # Infer link-mode:
  row <- unlist(lapply(X,nrow))
  col <- unlist(lapply(X,ncol))
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
  }
  mod$call <- match.call()
  class(mod) <- c('multiblock','list')
  return(mod)
}

#' @rdname unsupervised 
#' @export
gca <- function(X, ncomp=2, svd=TRUE, tol=10^-12, corrs=TRUE, ...){
  if(svd){
    obj <- gca.svd(X=X, tol=tol)
    obj$call = match.call()
    return(obj)
  } else {
    obj <- gca.rgcca(X=X, scale=FALSE, ncomp=ncomp, corrs=corrs, ...)
    obj$call = match.call()
    return(obj)
  }
}
gca.rgcca <- function(X, scale=FALSE, ncomp=1, corrs=TRUE, ...){
  n_block <- length(X)
  if(length(ncomp)==1){ ncomp <- rep(ncomp,n_block+1) }
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
gca.svd <- function(X, tol=10^-12){
  n_block <- length(X)
  n <- nrow(X[[1]])
  X <- lapply(X, function(i) scale(i, scale = FALSE))
  minrank <- min(min(unlist(lapply(X,ncol))),n)
  T <- Pb <- list()
  for(i in 1:n_block){
    udv <- svd(X[[i]])
    thisrank <- ifelse(sum(udv$d>tol) > 1, sum(udv$d>tol), 1)
    minrank <- min(thisrank, minrank)
    T[[i]]  <- udv$u[,1:thisrank,drop=FALSE]/rep(apply(udv$u[,1:thisrank,drop=FALSE],2,sd),each=n)
    Pb[[i]] <- udv$v[,1:thisrank,drop=FALSE] * rep(udv$d[1:thisrank], each=nrow(udv$v))
  }
  
  udv <- svd(do.call(cbind,T))
  C <- udv$u[,1:minrank,drop=FALSE]
  P <- udv$v[,1:minrank,drop=FALSE] * rep(udv$d[1:minrank], each=nrow(udv$v))
  A <- U <- list()
  for(i in 1:n_block){
    A[[i]] <- pracma::pinv(X[[i]]) %*% C
    U[[i]] <- X[[i]] %*% A[[i]]
  }
  
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
               blockScores = "Canonical scores", blockLoadings = "Canonical loadings")
  colnames(C) <- colnames(P) <- paste0('Comp ', 1:ncol(C))
  rownames(C) <- rownames(X[[1]])
  for(i in 1:length(X)){
    colnames(U[[i]])  <- paste0('Comp ', 1:ncol(U[[i]]))
    colnames(Pb[[i]]) <- paste0('Comp ', 1:ncol(Pb[[i]]))
    rownames(U[[i]])  <- rownames(X[[i]])
    rownames(Pb[[i]]) <- colnames(X[[i]])
  }
  rownames(P) <- unlist(lapply(Pb,rownames))
  obj <- list(scores=C, loadings=P, blockScores=U, blockLoadings=Pb, 
              blockCoef=A, cor=R, info=info)
  class(obj) <- list('multiblock','list')
  return(obj)
}

#' @rdname unsupervised 
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
              info = info, GPA = ret, call = match.call())
  class(obj) <- c('multiblock','list')
  return(obj)
}

#' @rdname unsupervised 
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
  attr(obj$scores, 'explvar') <- attr(obj$loadings, 'explvar') <- 100*(ret$global.pca$svd$vs^2/sum(ret$global.pca$svd$vs))
  for(i in 1:length(X)){
    explvar <- 100*(ret$separate.analyses[[i]]$svd$vs^2/sum(ret$separate.analyses[[i]]$svd$vs^2))
    attr(obj$blockScores[[i]], 'explvar') <- attr(obj$blockLoadings[[i]], 'explvar') <- explvar
    colnames(obj$blockScores[[i]]) <- colnames(obj$blockLoadings[[i]]) <- paste0('Comp ',1:ncol(obj$blockScores[[i]]))
  }
  class(obj) <- c('multiblock','list')
  return(obj)
}

#' @rdname unsupervised 
#' @export
pcagca <- function(X, commons=2, auto=TRUE, auto.par=list(explVarLim=40, rLim=0.8),
                   manual.par=list(ncomp=0, ncommon=0), tol=10^-12){
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
    U[[i]] <- udv$u[,1:ncomp[i],drop=FALSE]
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
              info=info, call=match.call())
  class(obj) <- c("multiblock","list")
  return(obj)
}

#' @rdname unsupervised 
#' @export
disco <- function(X, ncomp = 2, ...){
  ret <- RegularizedSCA::DISCOsca(do.call(cbind,X) , ncomp, unlist(lapply(X,ncol)))
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
              info = info, DISCOsca = ret, call = match.call())
  class(obj) <- c("multiblock","list")
  return(obj)
}

#' @rdname unsupervised 
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
  class(obj) <- c("multiblock","list")
  return(obj)
}

# Multiple Co-Inertia Analysis
#' @rdname unsupervised 
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
  return(list(A = res$astar, B = NULL, X = X, rgcca = res))
}

#' @rdname unsupervised
#' @export
jive <- function(X, ...){
  # X <- lapply(X,t) # Objects in columns
  r.jive::jive(X, ...)
  # jive(X, rankJ = 1, rankA = rep(1, length(data)), method = "perm",
  #      dnames = names(data), conv = "default", maxiter = 1000, scale = TRUE, center = TRUE,
  #      orthIndiv = TRUE, est = TRUE, showProgress=TRUE)
}

#' @rdname unsupervised 
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
  class(obj) <- c("multiblock","list")
  return(obj)
  # Has plot and print in ade4
  # Clustatis in https://cran.r-project.org/web/packages/ClustBlock/ClustBlock.pdf
}

#' @rdname unsupervised 
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
  class(obj) <- c("multiblock","list")
  return(obj)
}
