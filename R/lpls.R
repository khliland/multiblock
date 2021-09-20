#' L-PLS regression
#'
#' @param X1 \code{matrix} of size IxN (middle matrix)
#' @param X2 \code{matrix} of size IxJ (left matrix)
#' @param X3 \code{matrix} of size KxN (top matrix)
#' @param ncomp number of L-PLS components
#' @param doublecenter \code{logical} indicating if centering should be done both ways for X1 (default=TRUE)
#' @param scale \code{logical vector} of length three indicating if each of the matrices should be autoscaled.
#' @param type \code{character} indicating type of L-PLS ("exo"=default, "exo_ort" or "endo")
#' @param impute \code{logical} indicating if SVD-based imputation of missing data is required.
#' @param niter \code{numeric} giving number of iterations in component extraction loop.
#' @param subsetX2 \code{vector} defining optional sub-setting of X2 data.
#' @param subsetX3 \code{vector} defining optional sub-setting of X3 data. 
#' @param ... Additional arguments, not used.
#' 
#' @description Simultaneous decomposition of three blocks connected in an L pattern. 
#' 
#' @details Two versions of L-PLS 
#' are available: exo- and endo-L-PLS which assume an outward or inward relationship between the
#' main block X1 and the two other blocks X2 and X3.
#'
#' The \code{exo_ort} algorithm returns orthogonal scores and should be chosen for visual 
#' exploration in correlation loading plots. If exo-L-PLS with prediction is the main purpose 
#' of the model then the non-orthogonal \code{exo} type L-PLS should be chosen for which the 
#' predict function has prediction implemented.
#' 
#' ![](LPLSsmall.png "L-PLS diagram")
#'
#' @return An object of type \code{lpls} and \code{multiblock} containing all results from the L-PLS
#' analysis. The object type \code{lpls} is associated with functions for correlation loading plots, 
#' prediction and cross-validation. The type \code{multiblock} is associated with the default functions
#' for result presentation (\code{\link{multiblock_results}}) and plotting (\code{\link{multiblock_plots}}).
#' 
#' @author Solve Sæbø (adapted by Kristian Hovde Liland)
#' 
#' @references 
#' * Martens, H., Anderssen, E., Flatberg, A.,Gidskehaug, L.H., Høy, M., Westad, F.,Thybo, A., and Martens, M. (2005). Regression of a data matrix on descriptors of both its rows and of its columns via latent variables: L-PLSR. Computational Statistics & Data Analysis, 48(1), 103 – 123.
#' * Sæbø, S., Almøy, T., Flatberg, A., Aastveit, A.H., and Martens, H. (2008). LPLS-regression: a method for prediction and classification under the influence of background information on predictor variables. Chemometrics and Intelligent Laboratory Systems, 91, 121–132.
#' * Sæbø, S., Martens, M. and Martens H. (2010) Three-block data modeling by endo- and exo-LPLS regression. In Handbook of Partial Least Squares: Concepts, Methods and Applications. Esposito Vinzi, V.; Chin, W.W.; Henseler, J.; Wang, H. (Eds.). Springer.
#'
#' @examples
#' # Simulate data set
#' sim <- lplsData(I = 30, N = 20, J = 5, K = 6, ncomp = 2)
#' X1  <- sim$X1; X2 <- sim$X2; X3 <- sim$X3
#' lp  <- lpls(X1,X2,X3) # exo-L-PLS
#' 
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Functions for computation and extraction of results and plotting are found in \code{\link{lpls_results}}.
#' @export
lpls <- function(X1,X2,X3, ncomp = 2, doublecenter = TRUE, scale = c(FALSE,FALSE,FALSE),
                 type = c("exo"), impute = FALSE, niter = 25, subsetX2 = NULL, subsetX3 = NULL,...){
  
  # Rename inputs from Smilde, Liland and Næs 2021 to Sæbø et al.
  X1a <- X1
  X1 <- X2
  X2 <- X1a
  X3 <- t(X3)
  
  #The three matrices are assumed positioned in this manner
  #         Sæbø et al.               |       Smilde, Næs and Liland 2021
  #                 _________         |         _______ 
  #                |         |        |        |       |
  #                |         |        |        |       |
  #                |  t(X3)  |        |        |   X3  |
  #                |         |        |        |       |
  #                |_________|        |        |_______|
  #                                   |                   
  #                                   |                   
  #   _______       _________         |         _______       _________ 
  #  |       |     |         |        |        |       |     |         |
  #  |       |     |         |        |        |       |     |         |
  #  |  X1   |     |   X2    |        |        |  X1   |     |   X2    |
  #  |       |     |         |        |        |       |     |         |
  #  |_______|     |_________|        |        |_______|     |_________|
  #                           
  
  #  if(type!="exo_ort") impute=TRUE  
  if(impute){
    if(any(is.na(X1))) X1 <- svd.imp(X1)
    if(any(is.na(X2))) X2 <- svd.imp(X2)
    if(any(is.na(X3))) X3 <- svd.imp(X3)
  }
  
  # Extracting subsets
  if(!is.null(subsetX2)){
    X1 <- X1[subsetX2,,drop=F]
    X2 <- X2[subsetX2,,drop=F]
  }
  
  if(!is.null(subsetX3)){
    X3 <- X3[subsetX3,,drop=F]
    X2 <- X2[,subsetX3,drop=F]
  }
  
  # For saving 
  X1save <- X1
  X2save <- X2
  X3save <- X3
  
  # Dimensions
  X1dim <- dim(X1)
  X2dim <- dim(X2)
  X3dim <- dim(X3)
  
  # Centering and scaling
  X1  <- scale(X1, scale=scale[2])
  mX1 <- attr(X1, "scaled:center")
  X3  <- scale(X3, scale=scale[3])
  mX3 <- attr(X3, "scaled:center")
  
  rowmX2   <- apply(X2,1,mean, na.rm=TRUE)  
  colmX2   <- apply(X2,2,mean, na.rm=TRUE)
  grandmX2 <- mean(X2, na.rm=TRUE)    
  
  # Do double centering of X2?
  if(doublecenter){
    X2 <- X2 - t(matrix(rep(1, X2dim[2]), ncol = 1)%*%rowmX2) - 
      matrix(rep(1,X2dim[1]), ncol = 1)%*%colmX2 + 
      matrix(grandmX2, nrow=X2dim[1], ncol = X2dim[2])  
  } else {
    X2 <- X2 - matrix(grandmX2,nrow=X2dim[1],ncol=X2dim[2])
  }
  # Column scaling of X2
  X2 <- scale(X2, scale = scale[1])
  
  attributes(X1save) <- attributes(X1)
  attributes(X2save) <- attributes(X2)
  attributes(X3save) <- attributes(X3)  
  
  dlist <- list(X1 = as.matrix(t(X1)), X2 = as.matrix(X2), X3 = as.matrix(X3))
  
  
  T11 <- T12 <- T21 <- T22 <- T31 <- T32 <- numeric(0)
  W21 <- W22 <- numeric(0)
  P1  <- P3 <- P21 <- P22 <- numeric(0)
  D <- diag(rep(1,ncomp))
  X1totvar <- X2totvar <- X3totvar <- rep(0,ncomp+1)   
  
  for(a in 1:ncomp){
    X1totvar[a] <- sum(dlist$X1^2, na.rm=TRUE)
    X2totvar[a] <- sum(dlist$X2^2, na.rm=TRUE) 
    X3totvar[a] <- sum(dlist$X3^2, na.rm=TRUE)
    
    latent <- extractscores(dlist,niter=niter)
    
    # Score vectors   
    T11 <- cbind(T11, vecnorm(latent[[1]]$t1))
    T12 <- cbind(T12, vecnorm(latent[[1]]$t2))        
    T21 <- cbind(T21, vecnorm(latent[[2]]$t1))        
    T22 <- cbind(T22, vecnorm(latent[[2]]$t2))        
    T31 <- cbind(T31, vecnorm(latent[[3]]$t1))        
    T32 <- cbind(T32, vecnorm(latent[[3]]$t2))         
    
    
    #Fitting an exo-LPLS model
    if(type=="exo_ort"){
      #Orthogonal scores
      P1 <- cbind(P1, projectonto(dlist$X1,T21[,a]))
      P3 <- cbind(P3, projectonto(dlist$X3,T22[,a]))
      
      P21 <- cbind(P21, projectonto(dlist$X2,T21[,a]))
      P22 <- cbind(P22, projectonto(dlist$X2,T22[,a]))
      
      d <- drop(solve(crossprod(T21[,a]))%*%t(T21[,a])%*%P22[,a])
      D[a,a] <- d
      
      #Deflation
      dlist$X1 <- dlist$X1 - P1[,a]%*%t(T21[,a])
      dlist$X3 <- dlist$X3 - T22[,a]%*%t(P3[,a])  
      dlist$X2 <- dlist$X2 - T21[,a]%*%t(P21[,a]) - P22[,a]%*%t(T22[,a]) + T21[,a]%*%t(T22[,a])*drop(d)
    }#end type=="exo"
    
    
    if(type=="exo"){
      #Simpler model for X2, but non-orthogonal scores
      P1 <- projectonto(t(X1),T21)
      P3 <- projectonto(t(X3),T22)
      
      P21 <- projectonto(t(X2),T21)
      P22 <- projectonto(X2,T22)
      
      D <- solve(crossprod(T21))%*%t(T21)%*%P22
      
      #Deflation
      dlist$X1 <- t(X1 - T21%*%t(P1))
      dlist$X3 <- X3 - T22%*%t(P3) 
      dlist$X2 <- X2 - T21%*%D%*%t(T22)
      
    }#end type=="exo"        
    
    
    if(type=="endo"){
      P1 <- projectonto(t(X1),T12)
      P3 <- projectonto(t(X3),T31)
      D <- solve(crossprod(T12))%*%t(T12)%*%projectonto(X2, T31)
      
      #Deflation
      dlist$X1 <- t(X1 - T12%*%t(P1))
      dlist$X3 <- X3 - T31%*%t(P3)                         
      dlist$X2 <- X2 - T12%*%D%*%t(T31)
      
    }#end type=="endo"
    
  }
  
  X1totvar[ncomp+1] <-sum(dlist$X1^2, na.rm=TRUE)
  X2totvar[ncomp+1] <-sum(dlist$X2^2, na.rm=TRUE) 
  X3totvar[ncomp+1] <-sum(dlist$X3^2, na.rm=TRUE)
  
  X1varprop <- diff((X1totvar[1]-X1totvar)/(X1totvar[1]))
  X2varprop <- diff((X2totvar[1]-X2totvar)/(X2totvar[1]))
  X3varprop <- diff((X3totvar[1]-X3totvar)/(X3totvar[1]))
  
  
  
  #Various output
  B1 <- B3 <- Ca <- NULL
  if(type!="endo"){ 
    B1 <- T31%*%solve(t(P21)%*%T31)%*%t(P1)
    B3 <- T12%*%solve(t(P22)%*%T12)%*%t(P3)
#    options(warn=-1)
    suppressWarnings({
    R1  <- cor(X1,T21, use="pairwise.complete.obs") 
    R21 <- t(cor(T21,X2, use="pairwise.complete.obs")) 
    R22 <- t(cor(T22,t(X2), use="pairwise.complete.obs"))  
    R3  <- cor(X3,T22, use="pairwise.complete.obs")  
    R2rmean <- cor(rowmX2,T21, use="pairwise.complete.obs")  
    R2cmean <- cor(colmX2,T22, use="pairwise.complete.obs")
    })
#    options(warn=0)
  } else if(type=="endo"){
    V1 <- T11%*%solve(t(P1)%*%T11)
    V3 <- T32%*%solve(t(P3)%*%T32)
    Ca <- V1%*%D%*%t(V3)
#    options(warn=-1)
    suppressWarnings({
    R1  <- cor(X1,T12, use="pairwise.complete.obs") 
    R21 <- t(cor(T12,X2, use="pairwise.complete.obs")) 
    R22 <- t(cor(T31,t(X2), use="pairwise.complete.obs"))  
    R3  <- cor(X3,T31, use="pairwise.complete.obs")  
    R2rmean <- cor(rowmX2,T12, use="pairwise.complete.obs")  
    R2cmean <- cor(colmX2,T31, use="pairwise.complete.obs")
    })
#    options(warn=0)
  }
  
  res <- list(call=match.call())
  res$ncomp <- ncomp
  res$coefficients <- list(B1=B1, B3=B3, C=Ca)
  res$blockScores  <- list(T11=T11, T12=T12, T21=T21, T22=T22, T31=T31, T32=T32)
  res$blockLoadings<- list(P1=P1, P3=P3, P21=P21, P22=P22, D=D)
  res$corloadings  <- list(R1=R1, R21=R21, R22=R22, R3=R3, R2rmean=R2rmean, R2cmean=R2cmean)
  res$means     <- list(mX1=mX1, mX3=mX3, grandmX2=grandmX2, rowmX2=rowmX2, colmX2=colmX2)  
  res$data      <- list(X1=X1save, X2=X2save, X3=X3save)
  res$residuals <- dlist
  res$options   <- list(doublecenter=doublecenter, scale=scale, type=type)
  res$vars      <- list(X1varprop=X1varprop, X2varprop=X2varprop, X3varprop=X3varprop)
  res$info      <- list(method = "L-PLS", 
                        scores = "Not used", loadings = "Not used",
                        blockScores = "Block scores", blockLoadings = "Block loadings")
  class(res)    <- c('lpls', 'multiblock','list')
  return(res)
}



# ---------------------------------------------------------------
extractscores <- function(datalist, niter = 25, truncperc = NULL, truncvec = NULL){
  nmat  <- length(datalist)
  dims  <- lapply(datalist, dim)
  dimok <- TRUE
  for(ii in 1:(nmat-1)){
    dimok <- ifelse(dims[[ii]][2]==dims[[(ii+1)]][1], TRUE, FALSE)
  }
  if(!dimok) stop("Dimension mismatch\n")
  scorelist<-list()
  for(ii in 1:nmat){
    scorelist[[ii]]<-list(t1 = matrix(0, nrow = dims[[ii]][1], ncol = 1),
                          t2 = matrix(0, nrow = dims[[ii]][2], ncol = 1))
  }
  
  #Initiating the NIPALS algorithm
  scorelist[[1]]$t1 <- as.matrix(rnorm(dims[[1]][1], 0, 1))
  
  #NIPALS
  for(k in 1:niter){
    scorelist[[1]]$t2 <- projectonto(datalist[[1]], scorelist[[1]]$t1)
    for(j in 2:nmat){
      scorelist[[j]]$t2 <- projectonto(datalist[[j]], scorelist[[(j-1)]]$t2)
    }
    scorelist[[nmat]]$t1 <- projectonto(datalist[[nmat]], scorelist[[nmat]]$t2)
    for(j in (nmat-1):1){
      scorelist[[j]]$t1 <- projectonto(datalist[[j]], scorelist[[(j+1)]]$t1)
    }
  }
  return(scorelist)
}

# ------------------------------------------------------------

projectonto <- function(A,b){
  A.na <- any(is.na(A))  
  if(is.null(dim(A))) A <- matrix(A, ncol = 1)
  if(is.null(dim(b))) b <- matrix(b, ncol = 1)
  
  if(!any(dim(A)==dim(b)[1])) stop("Non-matching dimensions")
  q <- dim(b)[2]
  
  #Orthogonalize b if multiple b's
  if(A.na && q>1){
    svd1 <- svd(b)
    b <- b%*%svd1$v
  }
  
  if(dim(A)[1]==dim(b)[1]) A <- t(A)
  if(!A.na){ 
    proj <- A%*%b%*%solve(crossprod(b))
  } else {
    miss <- which(is.na(A))
    A[miss] <- 0    
    if(q==1){
      bb <- t(b^2)
      ones <- matrix(1,dim(A)[2],dim(A)[1])
      ones[miss] <- 0
      btbinv <- 1/drop(bb%*%ones)
      proj <- (A%*%b)*btbinv      
    } else if(q>1){
      proj <- numeric()
      for(i in 1:q){
        bb <- t(b[,i]^2)
        ones <- matrix(1,dim(A)[2],dim(A)[1])
        ones[miss] <- 0
        btbinv <- 1/drop(bb%*%ones)
        proj <- cbind(proj,(A%*%b[,i])*btbinv)      
      }
      proj <- proj%*%t(svd1$v)
    }
  }
  return(proj)
}


vecnorm <- function(vec){
  vec / sqrt(drop(crossprod(vec)))
}


svd.imp <- function(X, max.niter=50, expl.min=0.98, interactive=FALSE, tol=1e-3, ploteval=FALSE){
  
  ncomp = min(dim(X))
  missing <- which(is.na(X))
  X.scaled <- scale(X,scale=FALSE)
  colmeans <- attr(X.scaled,"scaled:center")
  meanmat <- matrix(1,nrow=dim(X)[1],ncol=1)%*%colmeans
  X.imp <- X
  impvals <- meanmat[missing]
  X.imp[missing] <- impvals
  initiate <- TRUE
  j<-1
  relchange <- 1
  change <- numeric()
  while(relchange > tol & j<=max.niter){        
    uvd <- svd(X.imp)
    
    if(initiate){
      ssx <- rep(0,ncomp)
      for(i in 1:ncomp){
        D <- matrix(0,i,i)
        diag(D) <- uvd$d[1:i]
        ssx[i] <- sum((uvd$u[,1:i,drop=F]%*%D%*%t(uvd$v[,1:i]))^2)
        initiate <- FALSE
      }
      ssxtot <- sum(X.imp^2)
      explvar <- ssx/ssxtot
      if(interactive){
        plot(1:ncomp,explvar)
        ncomp <- readline("Choose number of components for imputation \n")
        ncomp <- as.numeric(ncomp)
      }else{
        ncomp <- min(which(explvar > expl.min))  
      }
    }
    D <- matrix(0,ncomp,ncomp)
    diag(D) <- uvd$d[1:ncomp]        
    Xhat <- uvd$u[,1:ncomp,drop=F]%*%D%*%t(uvd$v[,1:ncomp,drop=F])
    impvallength <- sqrt(sum((impvals)^2))
    change[j] <- relchange <- sqrt(sum((impvals-Xhat[missing])^2))
    impvals <- Xhat[missing]
    X.imp[missing] <- impvals
    if(ploteval){
      plot(0:j,c(1,change/impvallength),ylab="Change",main=paste("Relative change in imputed values using",ncomp,"components\n"),xlab="iteration")
    }
    j <- j+1
  }
  X.imp
}
