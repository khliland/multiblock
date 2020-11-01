#' @name unsupervised
#' @title Unsupervised Multiblock Methods
#' @aliases hpca gca sca mfa mcoa pcagca disco statis
# # GSVD JIVE
#' @importFrom RGCCA rgcca
#' @importFrom MFAg MFA
#' @importFrom RegularizedSCA DISCOsca
#' @importFrom ade4 statis
#' 
#' @description Collection of unsupervised multiblock methods:
#' * Hierarchical Principal component analysis (\code{hpca})
#' * Generalized Canonical Analysis (\code{gca})
#' * Simultaneous Component Analysis (\code{sca})
#' * Higher Order Generalized SVD (\code{hogsvd})
#' * Multiple Factor Analysis (\code{mfa})
#' * Multiple Co-Inertia Analysis (\code{mcoa})
#' * PCA-GCA (\code{pcagca})
#' * Distinctive and Common Components with SCA (\code{disco})
#' * STATIS - Structuration des Tableaux A Trois Indices de la Statistique (\code{statis})
#' 
#' Original documentation of STATIS: \link[ade4]{statis}.
# Include also Penalized Exponential SCA (SLIDE, penalty-based)
# 
# NB! Removed:
#  r.jive jive

#' @rdname unsupervised 
#' @export
hpca <- function(X, ncomp=1, scale=FALSE, verbose=FALSE, ...){
  n_block <- length(X)
  if(length(ncomp)==1){ ncomp <- rep(ncomp,n_block) }
  X <- lapply(X, function(i) scale(i, scale = FALSE))
  X[[n_block+1]] <- do.call(cbind, X)
  C <- matrix(0, n_block+1, n_block+1)
  C[n_block+1,1:n_block] <- 1; C[1:n_block,n_block+1] <- 1
  res <- RGCCA::rgcca(A = X, C = C, tau=c(rep(1,n_block),0), scale = scale, ncomp=ncomp, scheme = function(x) x^4, verbose = verbose, ...)
  # A = scores and superscore (last), B = loadings and superloadings (last)
  return(list(A = lapply(1:(n_block+1), function(i)X[[i]]%*%res$astar[[i]]), B = res$astar, X = X, rgcca = res))
}

#' @rdname unsupervised 
#' @export
gca <- function(X, ncomp=1, svd=TRUE, tol=10^-12, corrs=TRUE, ...){
  if(svd){
    return(gca.svd(X=X, tol=tol))
  } else {
    gca.rgcca(X=X, scale=FALSE, ncomp=ncomp, corrs=corrs, ...)
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
  T <- list()
  for(i in 1:n_block){
    udv <- svd(X[[i]])
    thisrank <- ifelse(sum(udv$d>tol) > 1, sum(udv$d>tol), 1)
    minrank <- min(thisrank, minrank)
    T[[i]] <- udv$u[,1:thisrank,drop=FALSE]/rep(apply(udv$u[,1:thisrank,drop=FALSE],2,sd),each=n)
  }
  
  udv <- svd(do.call(cbind,T))
  C <- udv$u[,1:minrank,drop=FALSE]
  A <- U <- list()
  for(i in 1:n_block){
    A[[i]] <- pracma::pinv(X[[i]]) %*% C
    # A[[i]] <- tcrossprod(solve(crossprod(X[[i]])), X[[i]]) %*% C
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
  list(A=A, U=U, T=T, C=C, R=R)
}

#' @rdname unsupervised
#' @export
sca <- function(X, ncomp=1, scale=FALSE, samplelinked = 'auto', ...){
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
  } else {
    scores <- PCA$scores
    nvar   <- lapply(X, ncol)
    lookup <- unlist(lapply(1:length(X), function(i)rep(i, nvar[[i]])))
    loadings <- lapply(1:length(X), function(i)PCA$loadings[lookup==i,,drop=FALSE])
    names(loadings) <- names(X)
  }
  return(list(loadings=loadings, scores=scores))
}

#' @rdname unsupervised 
#' @export
mfa <- function(X, groupType = rep("n", length(X)), groupName = NULL, ...){
  MFAg::MFA(do.call(cbind,X), unlist(lapply(X,ncol)), TipoGrupo = groupType, NomeGrupos = groupName)
}

# Multiple Co-Inertia Analysis
#' @rdname unsupervised 
#' @export
mcoa <- function(X, ncomp=1, scale=FALSE, verbose=FALSE, ...){
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
pcagca <- function(X, commons=2, auto=TRUE, auto.par=list(explVarLim=40, rLim=0.8),
                   manual.par=list(ncomp=0, ncommon=0), tol=10^-12){
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
  
  # List of commons
  if(is.numeric(commons)){
    commons <- unique_combos(n_block, commons)
  }
  commonLabels <- lapply(commons, function(x)paste(x,sep="",collapse=","))
  
  # Basis scores
  for(i in 1:n_block){
    udv <- svd(X[[i]], ncomp[i], ncomp[i])
    U[[i]] <- udv$u[,1:ncomp[i],drop=FALSE]
    T[[i]] <- udv$u[,1:ncomp[i],drop=FALSE] %*% diag(udv$d[1:ncomp[i]])
    ncomp[i] <- pcaopt(X[[i]], T[[i]], udv$v, ncomp[i])
    T[[i]] <- T[[i]][,1:ncomp[i],drop=FALSE]
    blockLabels[[i]] <- character(0)
  }
  
  # Common components
  n_common <- length(commons)
  R <- numeric(0)
  for(i in 1:n_common){
    gcaComp <- gca(T[commons[[i]]])
    UC[[i]] <- gcaComp$U
    
    # Explained variance
    explVarX <- matrix(0, length(gcaComp$R),length(commons[[i]]))
    for(j in 1:length(commons[[i]])){
      XX <- X[[commons[[i]][j]]]
      for(k in 1:length(gcaComp$R)){
        ss <- gcaComp$U[[j]][,k]/c(sqrt(crossprod(gcaComp$U[[j]][,k])))
        loads <- crossprod(XX, ss)
        XX <- XX - tcrossprod(ss, loads)
        explVarX[k,j] <- (totVar[[commons[[i]][j]]]-sum(XX^2))/totVar[[commons[[i]][j]]]*100
      }
    }
    explVarX <- diff(rbind(numeric(length(commons[[i]])), explVarX))
    if(auto){
      nc <- which(gcaComp$R>auto.par$rLim & rowMeans(explVarX)>auto.par$explVarLim )
      ncommon[[i]] <- nc
    }

    if((!length(ncommon[[i]])==0) && ncommon[[i]]>0){
      R <- c(R, gcaComp$R[ncommon[[i]]])
      for(j in 1:length(commons[[i]])){
        # Store common scores
        u <- gcaComp$U[[j]][,ncommon[[i]],drop=FALSE]
        if(length(Scores)< commons[[i]][j]){
          Scores[[commons[[i]][j]]] <- matrix(0,nrow=n,ncol=0)}
        Scores[[commons[[i]][j]]] <- cbind(Scores[[commons[[i]][j]]], u/rep(sqrt(diag(crossprod(u))), each=n))
        for(k in 1:length(ncommon[[i]])){
          blockLabels[[commons[[i]][j]]] <- c(blockLabels[[commons[[i]][j]]], paste("C(", commonLabels[[i]],")-",k,sep=""))
        }
        # if nBlocks==2; Labels{idx(iblock)}=strcat('C_',num2str((1:length(nComps))')); else Labels{idx(iblock)}=strvcat(Labels{idx(iblock)},strcat(strcat('C', num2str(idx')','_'),num2str((1:length(nComps))'))); end

        # Deflate basis scores
        T[[commons[[i]][j]]] <- T[[commons[[i]][j]]] - u%*%solve(crossprod(u))%*%crossprod(u,T[[commons[[i]][j]]])
      }
    }
  }
  
  # Recompose basis into unique scores
  for(i in 1:n_block){
    udv <- svd(T[[i]])
    if(length(Scores)<i){
      Scores[[i]] <- matrix(0,nrow=n,ncol=0)}
    Scores[[i]] <- cbind(Scores[[i]], udv$u[,udv$d>tol, drop=FALSE])
    for(k in 1:sum(udv$d>tol)){
      blockLabels[[i]] <- c(blockLabels[[i]], paste("U-",k,sep=""))
    }
    # Labels{iblock}=strvcat(Labels{iblock},strcat('D',num2str(iblock),'_',num2str((1:dim_left)')));
  }
  
  # Calculate loadings
  for(i in 1:n_block){
    for(j in 1:ncol(Scores[[i]])){
      if(length(Loadings)<i){
        Loadings[[i]] <- matrix(0, nrow=ncol(X[[i]]), ncol=0)
      }
      Loadings[[i]] <- cbind(Loadings[[i]],crossprod(X[[i]], Scores[[i]][,j,drop=FALSE]))
      X[[i]] <- X[[i]] - tcrossprod(Scores[[i]][,j,drop=FALSE], Loadings[[i]][,j,drop=FALSE])
      if(length(explVar)<i){
        explVar[[i]] <- numeric(0)}
      explVar[[i]] <- c(explVar[[i]], (totVar[[i]]-sum(X[[i]]^2))/totVar[[i]]*100)
    }
  }
  names(ncommon) <- commonLabels
  explVar  <- lapply(1:n_block, function(i){e <- diff(c(0,explVar[[i]])); names(e) <- blockLabels[[i]];e})
  Scores   <- lapply(1:n_block, function(i){e <- Scores[[i]]; colnames(e) <- blockLabels[[i]];e})
  Loadings <- lapply(1:n_block, function(i){e <- Loadings[[i]]; colnames(e) <- blockLabels[[i]];e})
  if(auto){
    Scores   <- lapply(1:n_block,function(i)Scores[[i]][,explVar[[i]]>auto.par$explVarLim, drop=FALSE])
    Loadings <- lapply(1:n_block,function(i)Loadings[[i]][,explVar[[i]]>auto.par$explVarLim, drop=FALSE])
    explVar  <- lapply(1:n_block,function(i)explVar[[i]][explVar[[i]]>auto.par$explVarLim, drop=FALSE])
  }
  list(scores=Scores, loadings=Loadings, R=R, explVar=explVar, ncomp=unlist(lapply(explVar, length)), ncommon=ncommon)
}

#' @rdname unsupervised 
#' @export
disco <- function(X, ncomp, ...){
  RegularizedSCA::DISCOsca(do.call(cbind,X) , ncomp, unlist(lapply(X,ncol)))
}

# #' @rdname unsupervised 
# #' @export
# JIVE <- function(X, ...){
#   X <- lapply(X,t) # Objects in columns
#   r.jive::jive(X, ...)
#   # jive(X, rankJ = 1, rankA = rep(1, length(data)), method = "perm",
#   #      dnames = names(data), conv = "default", maxiter = 1000, scale = TRUE, center = TRUE,
#   #      orthIndiv = TRUE, est = TRUE, showProgress=TRUE)
# }

#' @rdname unsupervised 
#' @export
statis <- function(X, ncomp = 3, scannf = FALSE, tol = 1e-07, ...){
  K <- ktab.list.df(X)
  ade4::statis(K, scannf=scannf, nf=ncomp, tol=tol, ...)
  # Has plot and print in ade4
  # Clustatis in https://cran.r-project.org/web/packages/ClustBlock/ClustBlock.pdf
}

#' @rdname unsupervised 
#' @export
hogsvd  <- function(X, ncomp = 3){
  # TODO: Limit number of commponents in output
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
  
  # Eigen-decompose of S
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
  
  return(list(U=U, sigmas=sigmas, eigen_values=eigen_values, V=V, U_sigma=B))
}
