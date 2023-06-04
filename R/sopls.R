#' @title Sequential and Orthogonalized PLS (SO-PLS)
#' @importFrom Rcpp evalCpp
#' @importFrom pls plsr cvsegments
#' @importFrom SSBtools RowGroups
#' @importFrom progress progress_bar
#' @useDynLib multiblock
#' 
#' @param formula Model formula accepting a single response (block) and predictor block names separated by + signs.
#' @param ncomp Numeric vector of components per block or scalar of overall maximum components.
#' @param max_comps Maximum total number of components from all blocks combined (<= sum(ncomp)).
#' @param data The data set to analyse.
#' @param subset Expression for subsetting the data before modelling.
#' @param na.action How to handle NAs (no action implemented).
#' @param scale Logical indicating if variables should be scaled.
#' @param validation Optional cross-validation strategy "CV" or "LOO".
#' @param sequential Logical indicating if optimal components are chosen sequentially or globally.
#' @param segments Optional number of segments or list of segments for cross-validation. (See \code{[pls::cvsegments()]}).
#' @param sel.comp Character indicating if sequential optimal number of components should be chosen as minimum RMSECV ('opt', default) or by Chi-square test ('chi').
#' @param progress Logical indicating if a progress bar should be displayed while cross-validating.
#' @param ... Additional arguments to underlying methods.
#' 
#' @description Function for computing standard SO-PLS based on the interface of the \code{pls} package.
#' 
#' @details SO-PLS is a method which handles two or more input blocks by sequentially performing
#' PLS on blocks against a response and orthogonalising the remaining blocks on the extracted components.
#' Component number optimisation can either be done globally (best combination across blocks) or sequentially
#' (determine for one block, move to next, etc.).
#' 
#' @return An \code{sopls, mvr} object with scores, loadings, etc. 
#' associated with printing (\code{\link{sopls_results}}) and plotting methods (\code{\link{sopls_plots}}).
#' 
#' @references Jørgensen K, Mevik BH, Næs T. Combining designed experiments with several blocks of spectroscopic data. Chemometr Intell Lab Syst. 2007;88(2): 154–166.

#' @examples 
#' data(potato)
#' so <- sopls(Sensory ~ Chemical + Compression, data=potato, ncomp=c(10,10), 
#'             max_comps=10, validation="CV", segments=10)
#' summary(so)
#' 
#' # Scatter plot matrix with two first components from Chemical block
#' # and 1 component from the Compression block.
#' scoreplot(so, comps=list(1:2,1), ncomp=2, block=2)
#' 
#' # Result functions and more plots for SO-PLS 
#' # are found in ?sopls_results and ?sopls_plots.
#' @seealso SO-PLS result functions, \code{\link{sopls_results}}, SO-PLS plotting functions, \code{\link{sopls_plots}}, SO-PLS Måge plot, \code{\link{maage}}, and SO-PLS path-modelling, \code{\link{SO_TDI}}.
#' Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' @export
sopls <- function(formula, ncomp, max_comps = min(sum(ncomp), 20), data, 
                  subset, na.action, scale = FALSE, validation = c("none", "CV", "LOO"), 
                  sequential=FALSE, segments = 10, sel.comp='opt', progress=TRUE, ...){
  
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
  Y <- model.response(mf)
  response.type <- "continuous"
  if(is.factor(Y))
    response.type <- "categorical"
  if (is.matrix(Y)) {
    if (is.null(colnames(Y)))
      colnames(Y) <- paste("Y", 1:dim(Y)[2], sep = "")
  } else {
    if(is.factor(Y)){
      Y.dummy <- model.matrix(~Y-1, data.frame(Y = factor(Y)))
      Y <- as.matrix(as.numeric(levels(Y))[as.numeric(Y)])
    } else {
      Y <- as.matrix(Y)
      colnames(Y) <- deparse(formula[[2]])
    }
  }
  
  # # Convert factor blocks to dummy coding and make sure contents are matrices
  # M <- model.matrix(mt, mf)
  # blockColumns <- attr(M, 'assign')
  X <- lapply(X, function(x)if(is.factor(x)){return(dummycode(x))}else{return(x)})
  Xmeans <- Xscale <- list()
  for(i in 1:length(X)){
#   X[[i]] <- M[, blockColumns==i, drop=FALSE]
   Xmeans[[i]] <- apply(X[[i]],2,mean)
   if(scale){
     Xscale[[i]] <- apply(X[[i]],2,sd)
   }
  }
#  names(X) <- colnames(mf)[-1]
  X.concat <- do.call(cbind,X)
  # Check for missing dimnames
  if(is.null(rownames(X.concat)))
    rownames(X.concat) <- 1:nrow(X.concat)
  if(is.null(colnames(X.concat)))
    colnames(X.concat) <- 1:ncol(X.concat)
  
  y <- switch(response.type,
              continuous = Y,
              categorical = Y.dummy
  )
  
  ## Check components
  if(length(ncomp)==1 && length(X)>1){
    ncomp <- rep(ncomp, length(X))
  }
  ncomp <- pmin(ncomp, max_comps)
  
  ## Perform any scaling by sd per block:
  nobj <- dim(X[[1]])[1]
  
  ## Fit the SO-PLS model
  object <- sopls_modelling(X, Y, ncomp, max_comps, Xscale)
  object$fitted <- object$decomp$Ypred; object$decomp$Ypred <- NULL
  
  ## Validation
  Xs <- X
  if(scale){
    for(i in 1:length(Xs)){ # Centre, scale, de-centre (for cross-validation with segment-wise centring)
      Xs[[i]] <- ((Xs[[i]]-rep(Xmeans[[i]],each=nobj))/rep(Xscale[[i]],each=nobj))+rep(Xmeans[[i]],each=nobj)
    }
  }
  switch(match.arg(validation),
         CV = {
           if(sequential){
             val <- sopls_cv_seq(Xs, Y, ncomp, max_comps, sel.comp, segments, progress=TRUE, ...)
           } else {
             val <- sopls_cv(Xs, Y, ncomp, max_comps, segments, progress=TRUE, ...)
           }
         },
         LOO = {
           segments <- as.list(1:nobj)
           attr(segments, "type") <- "leave-one-out"
           if(sequential){
             val <- sopls_cv_seq(Xs, Y, ncomp, max_comps, sel.comp, segments, progress=TRUE, ...)
           } else {
             val <- sopls_cv(Xs, Y, ncomp, max_comps, segments, progress=TRUE, ...)
           }
         },
         none = {
           val <- NULL
         }
  )
  
  object$Xscale     <- Xscale
  object$na.action  <- attr(mf, "na.action")
  object$validation <- val
  object$call       <- match.call()
  object$model        <- mf
  object$Xmeans       <- Xmeans
  object$Ymeans       <- colMeans(Y)
  object$terms        <- mt
  object$ncomp        <- ncomp
  object$max_comps    <- max_comps
  class(object) <- c('sopls','mvr','multiblock')
  if(response.type == "categorical"){
    # object$classes <- sopls.classify(object, Y, X, ncomp, 'lda')
  }
  object
}

####################
# SO-PLS modelling #
####################
sopls_modelling <- function(X, Y, comps, max_comps = min(sum(comps), 20), Xscale){
  Y <- as.matrix(Y)
  n <- dim(Y)[1]
  nblock <- length(X)
  
  # Prepare XX''s
  C <- list()
  for(i in 1:nblock){
    X[[i]] <- as.matrix(X[[i]])
    if(length(Xscale)==0){
      C[[i]] <- tcrossprodQ(X[[i]] - rep(colMeans(X[[i]]), each=n))
    } else {
      C[[i]] <- tcrossprodQ((X[[i]] - rep(colMeans(X[[i]]), each=n))/rep(Xscale[[i]], each=n))
    }
  }
  
  so <- sopls_worker(C, Y, comps, max_comps, C, both=TRUE)
  list(decomp=so, data=list(X=X, Y=Y), method="PKPLS")
}


#####################
# SO-PLS prediction #
#####################
sopls_prediction <- function(SO, Xval, comps, scores){
  X <- SO$data$X
  Y <- SO$data$Y
  n    <- dim(X[[1]])[1]
  nval <- dim(Xval[[1]])[1]
  
  nblock  <- length(X)
  nresp   <- dim(Y)[2]
  selComp <- pathComp(comps, SO$decomp$compList)
  tot_comp <- length(selComp$hits)
  
  Y_pred <- array(0, c(nval, nresp, tot_comp))  
  Cr <- Crval <- 0
  for(i in 1:nblock){
    if(is.character(comps) || comps[i]>0){
      Xval[[i]] <- as.matrix(Xval[[i]])
      X[[i]]    <- as.matrix(X[[i]])
      Xval[[i]] <- Xval[[i]] - rep(SO$Xmeans[[i]], each = nval)
      X[[i]]    <- X[[i]]    - rep(SO$Xmeans[[i]], each = n)
      if(length(SO$Xscale)>0){
        Xval[[i]] <- Xval[[i]]/rep(SO$Xscale[[i]], each = nval)
        X[[i]]    <- X[[i]]   /rep(SO$Xscale[[i]], each = n)
      }
      Cr    <- Cr    + tcrossprodQ(X[[i]], X[[i]])    %*% SO$decomp$Ry[, selComp$hits, drop=FALSE]
      Crval <- Crval + tcrossprodQ(Xval[[i]], X[[i]]) %*% SO$decomp$Ry[, selComp$hits, drop=FALSE]
    }
  }
  # XW(P'W)^-1, ie. WB without Q
  no_Q <- Crval %*% solve(crossprodQ(SO$decomp$T[,selComp$hits,drop=FALSE], Cr))
  if(scores){
    dimnames(no_Q) <- list(rownames(X),apply(selComp$path,1,paste0, collapse=","))
    return(no_Q)
  }
  # Prediction per response
  for(r in 1:nresp){
    Yp_long <- t(apply(no_Q * rep(SO$decomp$Q[r,selComp$hits,drop=FALSE], each=nval), 1, cumsum))
    Y_pred[,r, ] <- Yp_long#[,(comp_curr-comp_last_block):comp_curr]
  }
  Y_pred <- Y_pred + rep(colMeans(Y), each=nval)
  dimnames(Y_pred) <- list(rownames(X),colnames(Y),apply(selComp$path,1,paste0, collapse=","))
  Y_pred
}


############################
# SO-PLS single prediction #  TODO: No scaling, consider trashing
############################
sopls_single_prediction <- function(X, Xval, Y, comps, max_comps = min(sum(comps), 20)){
  Y <- as.matrix(Y)
  n    <- dim(Y)[1]
  nval <- dim(Xval[[1]])[1]
  nblock <- length(X)
  
  # Prepare XX''s
  C <- Cval <- list()
  for(i in 1:nblock){
    Xval[[i]] <- as.matrix(Xval[[i]])
    X[[i]] <- as.matrix(X[[i]])
    Cval[[i]] <- tcrossprodQ(Xval[[i]] - rep(colMeans(X[[i]]), each=nval),X[[i]] - rep(colMeans(X[[i]]), each=n))
    C[[i]] <- tcrossprodQ(X[[i]] - rep(colMeans(X[[i]]), each=n))
  }
  
  sopls_worker(C, Y, comps, max_comps, Cval)
}

#########################
# Work-horse for SO-PLS #
#########################
# Performs either modelling/decomposition or
# efficient prediction, depending on inputs.
sopls_worker <- function(C, Y, comps, max_comps, Cval = NULL, both = FALSE){
  # Dimensions
  nblock <- length(C)
  n      <- dim(Y)[1]
  nresp  <- dim(Y)[2]
  
  # Create no-redundancy sequence of component
  cc <- component_combos(comps, max_comps)
  compList    <- cc$compList
  changeBlock <- cc$changeBlock
  blockIndex  <- cc$blockUsage$idx
  blockCombo  <- cc$blockUsage$groups
  nCombos     <- max(blockIndex)
  ncomps      <- rowSums(compList)
  tot_comps   <- length(changeBlock)
  
  # All combinations of block usage
  sumC <- list()
  for(i in 1:nCombos){
    sumC[[i]] <- 0
    for(j in 1:nblock){
      if(blockCombo[i,j]){
        sumC[[i]] <- sumC[[i]] + C[[j]]
      }
    }
  }
  
  # Check for prediction
  pred <- FALSE
  if(!is.null(Cval)){
    pred <- TRUE
    sumCval <- list()
    for(i in 1:nCombos){
      sumCval[[i]] <- 0
      for(j in 1:nblock){
        if(blockCombo[i,j]){
          sumCval[[i]] <- sumCval[[i]] + Cval[[j]]
        }
      }
    }
    Y_mean <- colMeans(Y)
    nval   <- dim(Cval[[1]])[1]
    Crval_currB <- list()
    for(b in 1:nblock){
      Crval_currB[[b]] <- matrix(0,nval,max_comps)
    }
  }
  
  # Prepare storage
  Ry      <- matrix(0.0, n, tot_comps)
  T       <- Ry
  Q       <- matrix(0.0, nresp, tot_comps)
  Ry_curr <- matrix(0.0, n, max_comps)
  T_curr  <- Ry_curr
  Q_curr  <- matrix(0.0, nresp, max_comps)
  Cr_currB<- list()
  Y_currB <- list()
  for(b in 1:nblock){
    Y_currB[[b]] <- Y
    Cr_currB[[b]] <- matrix(0,n,max_comps)
  }
  Y_curr  <- Y
  if(pred){
    Y_pred <- array(0.0, dim = c(nval,nresp,tot_comps))
  }
  
  # --------- Component extraction loop -------------
  for(comp in 2:tot_comps){
    cb <- changeBlock[comp]
    comp_curr <- ncomps[comp]
    Y_curr    <- Y_currB[[cb]]
    
    t <- C[[cb]] %*% Y_curr
    if(nresp > 1){
      w <- svd(crossprodQ(Y_curr,t))$v[,1]
      t <- t %*% w
    }
    if(comp_curr > 1){ # Orthogonalize on previous
      t <- t - T_curr[,1:(comp_curr-1),drop=FALSE] %*% crossprodQ(T_curr[,1:(comp_curr-1),drop=FALSE],t)
    }
    t <- t/sqrt(sum(t*t)) # FIXME: Save sum(t*t) for easy translation between ||t||=1 and ordinary t.
    if(nresp > 1){
      ry <- Y_curr %*% w
    } else {
      ry <- Y_curr
    }
    q <- crossprodQ(t,Y_curr)
    Y_curr <- Y_curr - t %*% q  # Deflation
    for(b in cb:nblock){
      Y_currB[[b]] <- Y_curr
    }
    
    # Store t, q, ry
    T_curr[, comp_curr]  <- t
    Q_curr[, comp_curr]  <- q
    Ry_curr[, comp_curr] <- ry
    if(both || !pred){
      T[, comp]  <- t
      Q[, comp]  <- q
      Ry[, comp] <- ry
    }
    if(pred){
      # Update "X_val*W" ~= C*Ry with current component
      Cr_currB[[cb]][,comp_curr]    <- C[[cb]] %*% ry
      Crval_currB[[cb]][,comp_curr] <- Cval[[cb]] %*% ry 
      if(cb < nblock){
        for(b in (cb+1):nblock){
          Cr_currB[[b]]    <- Cr_currB[[cb]]
          Crval_currB[[b]] <- Crval_currB[[cb]]
        }
      }
    }
    
    # Perform prediction at the end of each "curr"-series
    if(pred && (comp == tot_comps || changeBlock[comp+1] < nblock)){
      comp_last_block <- compList[comp,nblock] # Length of current series
      if(comp_curr-comp_last_block == 0){      # Compensate for first series starting at 0
        comp_last_block <- comp_last_block-1
      }
      
      # XW(P'W)^-1, ie. WB without Q
      no_Q <- Crval_currB[[cb]][,1:comp_curr,drop=FALSE] %*% 
        solve(crossprodQ(T_curr[,1:comp_curr,drop=FALSE], Cr_currB[[cb]][,1:comp_curr,drop=FALSE]))
      # Prediction per response
      for(r in 1:nresp){
        if(comp_curr==1){
          Yp_long <- no_Q * rep(Q_curr[r,1:comp_curr,drop=FALSE], each=nval)
        } else {
          Yp_long <- t(apply(no_Q * rep(Q_curr[r,1:comp_curr,drop=FALSE], each=nval), 1, cumsum))
        }
        Y_pred[,r, (comp-comp_last_block):comp] <- Yp_long[,(comp_curr-comp_last_block):comp_curr]
      }
    }
  } # -------------- End component extraction loop ---------------
  
  # If prediction, return predicted values
  if(both){
    Y_pred <- Y_pred + rep(Y_mean,each=nval)
    dimnames(Y_pred) <- list(rownames(Cval[[1]]), colnames(Y), apply(compList,1,paste, collapse = ","))
    cols <- apply(compList,1,paste, collapse = ",")
    dimnames(T) <- list(rownames(Y), cols)
    dimnames(Ry) <- list(rownames(Y), cols)
    dimnames(Q) <- list(colnames(Y), cols)
    return(list(Ypred=Y_pred, Q=Q, T=T, Ry=Ry, compList=compList, changeBlock=changeBlock))
  } else {
    if(pred){
      Y_pred <- Y_pred + rep(Y_mean,each=nval)
      dimnames(Y_pred) <- list(rownames(Cval[[1]]), colnames(Y), apply(compList,1,paste, collapse = ","))
      return(list(Ypred=Y_pred, compList=compList, changeBlock=changeBlock))
    } else {
      
      # Otherwise, return decomposition
      # TODO: Include W and P instead of Ry?
      cols <- apply(compList,1,paste, collapse = ",")
      dimnames(T) <- list(rownames(Y), cols)
      dimnames(Ry) <- list(rownames(Y), cols)
      dimnames(Q) <- list(colnames(Y), cols)
      return(list(Q=Q, T=T, Ry=Ry, compList=compList, changeBlock=changeBlock))
    }
  }
}

#############################################
# Create no-redundance sequence of component
component_combos <- function(comps, max_comps){
  nblock <- length(comps)
  # Determine block order
  unfiltered        <- rev(expand.grid(lapply(rev(comps), function(i) 0:i)))
  names(unfiltered) <- paste0('block ', 1:nblock)
  filtered    <- as.matrix(unfiltered[rowSums(unfiltered) <= max_comps,,drop=FALSE])
  changeBlock <- apply(diff(filtered),1,function(i)which(i>0)[1])
  
  # Determine involved blocks
  rg <- RowGroups(filtered!=0, TRUE)
  list(compList = filtered, changeBlock = unname(c(nblock,changeBlock)), blockUsage = rg)
}


#####################
# Path through compList
pathComp <- function(comps, compList){
  if(length(comps)==1 && comps=="all"){
    nb <- ncol(compList)
    ord <- order(as.numeric(apply(compList[, nb:1, drop=FALSE],1,paste, collapse="")))[-1]
    return(list(path = compList[ord,,drop=FALSE], hits = ord))
  }
  if(length(comps)==1 && comps=="raw"){
    ord <- 2:nrow(compList)
    return(list(path = compList[ord,,drop=FALSE], hits = ord))
  }
  nblocks <- length(comps)
  mat <- matrix(0, 0, nblocks)
  for(b in 1:nblocks){
    if(comps[b]>0){
      base <- 1:comps[b]
      if(b>1){
        for(c in seq(b-1,1)){
          base <- cbind(comps[c],base)
        }
      }
      if(b<nblocks){
        for(c in (b+1):nblocks){
          base <- cbind(base,0)
        }
      }
      mat <- rbind(mat,base)
    }
  }
  list(path = mat, hits = match(apply(mat,1,paste0, collapse=","), apply(compList,1,paste0, collapse=",")))
}


##############
# Chi2 error
chi2cv <- function(RMSECV, n, pvalue = 0.05){
  cvsig <- pchisq((n*min(RMSECV^2))/(RMSECV^2), n)
  ind <- which(cvsig > pvalue)[1]
  if(is.na(ind)){
    return(unname(which.min(RMSECV)))
  } else {
    return(unname(ind))
  }
}


###################
# Function testing
# st1 <- system.time(som <- sopls_modelling(list(halvdikotom$chemical, halvdikotom$NIR), halvdikotom$salt, c(2,10), max_comps = 10))
# st2 <- system.time(sosp <- sopls_single_prediction(list(halvdikotom$chemical[1:30,], halvdikotom$NIR[1:30,]), 
#                                             list(halvdikotom$chemical[-(1:30),], halvdikotom$NIR[-(1:30),]), 
#                                             halvdikotom$salt[1:30,], c(2,5), max_comps = 10))
# st3 <- system.time(cvp <- sopls_cv(list(halvdikotom$chemical, halvdikotom$NIR), halvdikotom$salt, c(2,5), 4, validation = "CV", k=10, type="consecutive"))
# 
# for(i in 1:5){
#   process_data[[i]] <- scale(process_data[[i]])
# }
# st1 <- system.time(som  <- sopls_modelling(process_data[c(1,4)], process_data[[5]], c(4,3), max_comps = 10))
# st1_f <- system.time(sos1f <- sopls_modelling(lapply(process_data[c(1,4)], function(i)i[-1,]),
#                                                   process_data[[5]][-1,], c(4,3), max_comps = 10))
# st2 <- system.time(sosp <- sopls_single_prediction(process_data[c(1,4)], process_data[c(1,4)], process_data[[5]], c(4,3), max_comps = 10))
# st_f <- system.time(sosf <- sopls_single_prediction(lapply(process_data[c(1,4)], function(i)i[-1,]), lapply(process_data[c(1,4)], function(i)i[1,,drop=FALSE]),
#                                                     process_data[[5]][-1,], c(4,3), max_comps = 10))
# st3 <- system.time(cvp <- sopls_cv(process_data[c(1,4)], process_data[[5]], c(4,3), 10, validation = "LOO", type="consecutive"))
# 
# st3CV <- system.time(cvpCV <- sopls_cv(process_data[c(1,4)], process_data[[5]], c(4,3), 10, validation = "CV", type="consecutive", k=10))
# 
# stseq <- sopls_cv_seq(process_data[c(1,4)], process_data[[5]], c(4,3), 10, validation = "CV", type="consecutive", k=10, sel.comp = "chi")
# 
# 
# ######################################
# X <- cbind(process_data[[1]], process_data[[4]]); Xt <- X[1,,drop=FALSE]; X <- X[-1,]
# Xc <- X-rep(colMeans(X),each=794)
# Xtc <- Xt-colMeans(X)
# Ry <- sos1f$decomp$Ry
# T <- sos1f$decomp$T
# W <- crossprod(Xc,Ry)
# for(i in 1:20){
#   if(sos1f$decomp$changeBlock[i]==1){
#     W[7:10,i] <- 0
#   } else {
#     W[1:6,i] <- 0
#   }
# }
# for(i in 1:20)W[,i] <- W[,i]/c(sqrt(crossprod(W[,i])))
# P <- crossprod(Xc,T)
# 
# pc <- pathComp(c(3,2), cvp$compList)
# B <- W[,pc$hits]%*%solve(t(P[,pc$hits])%*%W[,pc$hits])%*%t(Q[,pc$hits])
# 
# pred <- Xtc%*%B
# pred + colMeans(process_data[[5]][-1,])
