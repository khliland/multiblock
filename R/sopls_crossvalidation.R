######################################
# SO-PLS sequential cross-validation #
######################################
# SO-PLS sequential cross-validation

sopls_cv_seq <- function(X, Y, comps, max_comps, sel.comp, segments, progress, ...){
  Y <- as.matrix(Y)
  n      <- dim(Y)[1]
  nresp  <- dim(Y)[2]
  nblock <- length(X)
  
  socv <- sopls_cv(X=X, Y=Y, comps=comps, max_comps=max_comps, segments=segments, progress=progress, ...)
  RMSECV <- socv$RMSECV
  compList <- socv$compList
  rownames(compList) <- names(RMSECV)
  chosen <- numeric(nblock)
  
  # Automatic optimisation of number of components
  if(is.character(sel.comp)){
    for(b in 1:nblock){
      if(sel.comp == "opt"){
        if(b < nblock){
          chosen[b] <- which.min(RMSECV[rowSums(compList[,-1,drop=FALSE])==0])-1
          RMSECV   <- RMSECV[compList[,1]==chosen[b]]
          compList <- compList[compList[,1]==chosen[b],-1]
        } else {
          chosen[b] <- which.min(RMSECV)-1
        }
      } else {
        if(b < nblock){
          chosen[b] <- chi2cv(RMSECV[rowSums(compList[,-1,drop=FALSE])==0], n, 0.05)-1
          RMSECV   <- RMSECV[compList[,1]==chosen[b]]
          compList <- compList[compList[,1]==chosen[b],-1]
        } else {
          chosen[b] <- which.min(RMSECV)-1
        }
      }  
    }
    socv$chosen <- chosen
  } else {
    socv$chosen <- sel.comp
  }
  socv
}


###########################
# SO-PLS cross-validation #
###########################
sopls_cv <- function(X, Y, comps, max_comps, segments, progress=TRUE, ...){
  Y <- as.matrix(Y)
  n      <- dim(Y)[1]
  nresp  <- dim(Y)[2]
  nblock <- length(X)
  
  # Convert to correspondance vector
  mf <- match.call(expand.dots = TRUE)
  if(length(segments)==1 || !is.null(mf$segment.type)){
    segments <- cvsegments(n, k = segments, type = ifelse(is.null(mf$segment.type), 'random', mf$segment.type))
  }
  cv <- unlist(lapply(1:length(segments), function(i)rep(i,length(segments[[i]]))))
  nseg <- max(cv)

  # Initialize progress indicator
  if(progress){
    pb <- progress_bar$new(format = "  SO-PLS :what [:bar] :percent eta: :eta",
                           clear = TRUE, total = nseg)
    nwhite <- nchar(paste0("segment ", nseg))-nchar(paste0("nblock"))
    nnseg  <- nchar(paste0(nseg))
  }

  # Prepare XX''s end other CV accessories
  C <- Xc <- sumX <- mm <- Xm <- list()
  for(i in 1:nblock){
    # if(progress){
    #   pb$tick(tokens = list(what = paste0("block ", i)))
    # }
    X[[i]] <- as.matrix(X[[i]])
    if(length(segments) == n){
      Xc[[i]]  <- (rep(colSums(X[[i]]),each=n)-X[[i]])/(n-1) # Means without i-th sample
      Yc       <- (rep(colSums(Y),each=n)-Y)/(n-1)           # Means without i-th sample
    } else {
      segLength <- numeric(nseg)
      cvMat <- matrix(0,nseg,n)
      for(s in 1:nseg){
        segLength[s] <- sum(cv==s)
        cvMat[s,cv==s] = 1
      }
      sumX[[i]] <- colSums(X[[i]])
      sumy <- colSums(Y)
      Xc[[i]] <- (rep(sumX[[i]],each=nseg) - cvMat %*% X[[i]])/(n-segLength) # Means without i-th sample set
      Yc      <- (rep(sumy,each=nseg) - cvMat %*% Y)/(n-segLength)           # Means without i-th sample set
    }
    mm[[i]]  <- rowSums(Xc[[i]]^2)          # Inner products of Xc
    Xm[[i]]  <- tcrossprodQ(X[[i]],Xc[[i]]) # X*mean(X) for each sample mean
    C[[i]] <- tcrossprodQ(X[[i]])
  }
  
  n2 <- rep(FALSE, n)
  for(i in 1:nseg){
    if(progress){
      pb$tick(tokens = list(what = paste0("segment ", format(i, width=nnseg, justify="right"))))
    }
    inds2 <- n2
    inds2[cv==i] <- TRUE
    # Compute X*X' with centred X matrices excluding observation i
    Ci <- Cival <- Yi <- list()
    for(blk in 1:nblock){
      Ci[[blk]] <- C[[blk]] - (Xm[[blk]][,i] + matrix(Xm[[blk]][,i,drop=FALSE] - mm[[blk]][i], n,n, byrow=TRUE))
      Ci[[blk]][inds2,] <- 0
      nt <- sum(cv==i)
      # Compute Xval*Xval' with centred X matrices excluding observation i
      XXvt <- C[[blk]][cv==i,cv!=i] - (matrix(Xm[[blk]][cv!=i,i],nt,n-nt,byrow=TRUE) + (Xm[[blk]][i,i] - mm[[blk]][i]))
      Cival[[blk]] <- matrix(0,nt,n)
      Cival[[blk]][,!inds2] <- XXvt
      Ci[[blk]][,inds2] <-  0
    }
    Yi <- Y-rep(Yc[i,,drop=FALSE], each=n)
    Yi[inds2,] <- 0
    
    # Prediction
    cvresi <- sopls_worker(Ci, Yi, comps, max_comps, Cival)
    if(i == 1){
      Y_pred <- array(0.0, dim=c(n,nresp,dim(cvresi$Ypred)[3]))
    }
    Y_pred[inds2,,] <- cvresi$Ypred + rep(Yc[i,],each=nt)
  }
  dimnames(Y_pred) <- list(rownames(Y), colnames(Y), dimnames(cvresi$Ypred)[[3]])
  
  if(progress){
    pb$terminate()
  }
  RMSECV_ind <- sqrt(apply((Y_pred-rep(Y,dim(Y_pred)[3]))^2,2:3,mean))
  dimnames(RMSECV_ind) <- dimnames(Y_pred)[2:3]
  RMSECV     <- sqrt(apply(RMSECV_ind^2,2,mean))
  # vars       <- apply(Y,2,var)*(n-1)/n
  msep0 <- apply(Y,2,var)*(n-1)/n
  expl_var_ind <- 1 - RMSECV_ind^2/msep0 # RMSECV_ind[,1]^2 # vars
  expl_var     <- 1 - RMSECV^2 / mean(msep0) # RMSECV[1]^2 # mean(vars)
  list(Ypred=Y_pred, expl_var=expl_var, RMSECV=RMSECV, expl_var_ind=expl_var_ind, RMSECV_ind=RMSECV_ind, compList=cvresi[[2]], segments=segments)
}




