#' @name rosa
#' @title Response Oriented Sequential Alternation - ROSA
#'
#' @param formula Model formula accepting a single response (block) and predictor block names separated by + signs.
#' @param ncomp The maximum number of ROSA components.
#' @param Y.add Optional response(s) available in the data set.
#' @param common.comp Automatically create all combinations of common components up to length \code{common.comp} (default = 1).
#' @param data The data set to analyse.
#' @param subset Expression for subsetting the data before modelling.
#' @param na.action How to handle NAs (no action implemented).
#' @param scale Optionally scale predictor variables by their individual standard deviations.
#' @param weights Optional object weights.
#' @param validation Optional cross-validation strategy "CV" or "LOO".
#' @param internal.validation Optional cross-validation for block selection process, "LOO", "CV3", "CV5", "CV10" (CV-number of segments), or vector of integers (default = FALSE).
#' @param fixed.block integer vector with block numbers for each component (0 = not fixed) or list of length <= ncomp (element length 0 = not fixed).
#' @param design.block integer vector containing block numbers of design blocks
#' @param canonical logical indicating if canonical correlation should be use when calculating loading weights (default), enabling B/W maximization, common components, etc. Alternatively (FALSE) a PLS2 strategy, e.g. for spectra response, is used.
#' @param ... Additional arguments for \code{cvseg} or \code{rosa.fit}
#'
#' @return An object of classes \code{rosa} and \code{mvr} having several associated printing (\code{\link{rosa_results}}) and plotting methods (\code{\link{rosa_plots}}).
#' 
#' @description Formula based interface to the ROSA algorithm following the style of the \code{pls} package.
#' 
#' @details ROSA is an opportunistic method sequentially selecting components from whichever block explains
#' the response most effectively. It can be formulated as a PLS model on concatenated input block with block
#' selection per component. This implementation adds several options that are not described in the literature.
#' Most importantly, it opens for internal validation in the block selection process, making this more 
#' robust. In addition it handles design blocks explicitly, enables classification and secondary responses (CPLS),
#' and definition of common components.
#' 
#' @importFrom pracma mrdivide Rank
#' @importFrom RSpectra svds
#' @importFrom pls coefplot corrplot cvsegments loadingplot loading.weights loadings predplot scores scoreplot RMSEP validationplot R2 MSEP mvrValstats
#' @importFrom grDevices rgb
#' @importFrom graphics axis box lines par points text
#' @importFrom stats approx coef cor drop.terms fitted formula median model.frame model.matrix model.response na.pass napredict optimize pchisq predict rnorm sd terms update var
#' 
#' @export coefplot corrplot cvsegments loadingplot loading.weights loadings predplot scores scoreplot RMSEP validationplot R2 MSEP mvrValstats
#'
#' @references Liland, K.H., Næs, T., and Indahl, U.G. (2016). ROSA - a fast extension of partial least squares regression for multiblock data analysis. Journal of Chemometrics, 30, 651–662, doi:10.1002/cem.2824.
#' @examples
#' data(potato)
#' mod <- rosa(Sensory[,1] ~ ., data = potato, ncomp = 10, validation = "CV", segments = 5)
#' summary(mod)
#' 
#' # For examples of ROSA results and plotting see 
#' # ?rosa_results and ?rosa_plots.
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results and plotting are found in \code{\link{rosa_results}} and \code{\link{rosa_plots}}, respectively.
#' @export
rosa <- function(formula, ncomp, Y.add, common.comp = 1, data,
                 subset, na.action, scale = FALSE, weights = NULL,
                 validation = c("none", "CV", "LOO"), internal.validation = FALSE, fixed.block = NULL,
                 design.block = NULL, canonical = TRUE, ...){
  # Not implemented:
  # @param lower numerical indicating lower bound for powered PLS. Default = 0.5.
  # @param upper numerical indicating upper bound for powered PLS. Default = 0.5.
  
  ## Warn if impossible combinations
  if(!canonical && common.comp > 1){
    common.comp <- 1
    warning('Common components are not supported when PLS2 strategy is used (canonical=FALSE).')
  }
  ## Get the model frame
  mf <- match.call(expand.dots = FALSE)
  if (!missing(Y.add)) {
    ## Temporarily add Y.add to the formula
    Y.addname <- as.character(substitute(Y.add))
    mf$formula <- update(formula, paste("~ . +", Y.addname))
  }

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

  if (missing(Y.add)) {
    Y.add <- NULL
  } else {
    Y.add <- mf[,Y.addname]
    ## Remove Y.add from the formula again
    mt <- drop.terms(mt, which(attr(mt, "term.labels") == Y.addname),
                     keep.response = TRUE)
  }
  # X.concat <- pls:::delete.intercept(model.matrix(mt, mf))
  # X <- model.frame(delete.response(mt),mf)
  # Convert factor blocks to dummy coding and make sure contents are matrices
#  M <- model.matrix(mt, mf)
#  blockColumns <- attr(M, 'assign')
#  X <- list()
#  for(i in 1:(length(mf)-1)){
#    X[[i]] <- M[, blockColumns==i, drop=FALSE]
#  }
  # for(i in 1:length(X)){
  #   if(is.factor(X[[i]])){
  #     X[[i]] <- I(dummycode(X[[i]]))
  #   }
  #   X[[i]] <- I(as.matrix(X[[i]]))
  # }
  # names(X) <- colnames(mf)[-1]
  X <- lapply(X, function(x)if(is.factor(x)){return(dummycode(x))}else{return(x)})
  X.concat <- do.call(cbind,X)

  y <- switch(response.type,
              continuous = Y,
              categorical = Y.dummy
  )

  # Check ranks vs ncomp
  cRank <- 0
  for(i in 1:length(X)){
    cRank <- cRank + Rank(X[[i]])
  }
  nobj <- dim(X.concat)[1]
  if(cRank < ncomp || ncomp > (nobj-1)){
    ncomp <- min(cRank, nobj-1)
    warning(paste("ncomp reduced to ", ncomp, " because of predictor rank", sep=""))
  }

  ## Perform any scaling by sd per block:
  nobj <- dim(X[[1]])[1]
  if (scale) {
    for(i in 1:length(X)){
      # nobj <- dim(X[[i]])[2]
      ## This is faster than sd(X), but cannot handle missing values:
      the_scale <- sqrt(colSums((X[[i]] - rep(colMeans(X[[i]]), each = nobj))^2) /
                      (nobj - 1))
      if (any(abs(the_scale) < .Machine$double.eps^0.5))
        warning("Scaling with (near) zero standard deviation")
      X[[i]] <- X[[i]] / rep(the_scale, each = nobj)
    }
    # nobj <- dim(X.concat)[2]
    ## This is faster than sd(X), but cannot handle missing values:
    the_scale <- sqrt(colSums((X.concat - rep(colMeans(X.concat), each = nobj))^2) /
                    (nobj - 1))
    if (any(abs(the_scale) < .Machine$double.eps^0.5))
      warning("Scaling with (near) zero standard deviation")
    X.concat <- X.concat / rep(the_scale, each = nobj)
  }

  ## Convert fixed.block to list and check if all components are fixed
  fixed.order <- NULL
  if(!is.null(fixed.block)){
    if(is.numeric(fixed.block)){
      fb <- lapply(lapply(fixed.block, function(i)i), function(i)ifelse(i==0, return(NULL),return(i)))
      if(length(fixed.block) == ncomp && sum(fixed.block == 0) == 0){ # All fixed
        fixed.order <- fb
      }
      fixed.block <- fb
    } else {
      if(length(fixed.block) == ncomp && !any(unlist(lapply(fixed.block,is.null)))){ # All fixed
        fixed.order <- fixed.block
      }
    }
  }
  
  # Prepare internal validation for component order decision
  if(is.null(internal.validation) || isFALSE(internal.validation)){
    internal.validation <- FALSE
  } else {
    if(isTRUE(internal.validation) || any(internal.validation == "LOO"))
      internal.validation <- 1:nobj
    else {
      if(any(c("CV3","CV5","CV10")%in%internal.validation)){
        if(internal.validation == "CV3"){
          internal.validation <- rep(1:3,each=ceiling(nobj/3))[1:nobj]
        } else {
          if(internal.validation == "CV5"){
            internal.validation <- rep(1:5,each=ceiling(nobj/5))[1:nobj]
          } else {
            if(internal.validation == "CV10"){
              internal.validation <- rep(1:10,each=ceiling(nobj/10))[1:nobj]
            }
          }
        }
      } else {
        if(length(internal.validation) != nobj)
          stop('Internal validaiton vector must match number of objects.')
      }
    }
  }
  

  ## Fit the ROSA model
  object <- rosa.fit(X, X.concat, y, Y.add, ncomp, common.comp, weights, fixed.order, fixed.block, design.block, canonical, internal.validation, ...)

  ## Validation
  switch(match.arg(validation),
         CV = {
           val <- rosaCV(X, y, Y, ncomp, Y.add = Y.add, response.type = response.type, common.comp,
                         scale = scale, weights = weights, fixed.order = object$order, canonical = canonical, 
                         internal.validation = FALSE, ...)
         },
         LOO = {
           segments <- as.list(1:nobj)
           attr(segments, "type") <- "leave-one-out"
           val <- rosaCV(X, y, Y, ncomp, response.type = response.type, Y.add = Y.add, common.comp,
                         scale = scale, weights = weights, segments = segments, fixed.order = object$order,
                         canonical = canonical, internal.validation = FALSE, ...)
         },
         none = {
           val <- NULL
         }
  )

  if (is.numeric(scale)) object$scale <- scale
  object$na.action  <- attr(mf, "na.action")
  object$validation <- val
  object$call       <- match.call()
  object$design.block <- design.block
  object$canonical    <- canonical
  object$fixed.block  <- fixed.block
  object$model        <- mf
  object$X            <- X
  object$terms        <- mt
  class(object) <- c('rosa','mvr','multiblock')
  if(response.type == "categorical"){
    object$classes <- rosa.classify(object, Y, X, ncomp, 'lda')
  }
  object
}


#########################
# Internal ROSA fitting # ---------------------------------------------------------
#########################
rosa.fit <- function(X, X.concat, y, Y.add, ncomp, common.comp, weights, fixed.order, fixed.block, design.block, canonical, internal.validation, lower = 0.5, upper = 0.5, ...){
  nobj   <- dim(X.concat)[1]
  npred  <- dim(X.concat)[2]
  nresp  <- dim(y)[2]
  nadd   <- ifelse(is.null(Y.add),0,dim(Y.add)[2])
  nblock <- length(X)

  # Block sizes and centering
  m <- integer(nblock+1)
  if(is.null(weights)){
    for(i in 1:nblock){
      if(is.vector(X[[i]]) || is.null(dim(X[[i]]))){
        m[i] <- 1
        X[[i]] <- X[[i]] - rep(mean(X[[i]]), each = nobj)
      } else {
        m[i]   <- dim(X[[i]])[2]
        X[[i]] <- X[[i]] - rep(colMeans(X[[i]]), each = nobj)
      }
    }
    Xmeans   <- colMeans(X.concat)
  } else {
    for(i in 1:nblock){
      if(is.vector(X[[i]])){
        m[i] <- 1
      } else {
        m[i]   <- dim(X[[i]])[2]
      }
      X[[i]] <- X[[i]] - rep(crossprod(weights,X[[i]])/sum(weights), each = nobj)
    }
    Xmeans   <- crossprod(weights,X.concat)/sum(weights)
  }
  X.orig   <- X.concat
  X.concat <- X.concat - rep(Xmeans, each = nobj)
  Ymeans <- colMeans(y)
  y.orig <- y
  y <- y-rep(Ymeans, nobj)
  m[nblock+1] <- sum(m[-(nblock+1)])

  # Block combinations
  combos <- common.comb(nblock, common.comp)
  ncomb  <- nrow(combos)

  order <- list()            # Block order
  count <- integer(nblock+1) # Total block usage

  # Design block handling, maximum k-1 component per block
  if(!is.null(design.block)){
    ncomp.design <- integer(length(design.block))
    for(i in 1:length(design.block)){
      ncomp.design[i] <- Rank(X[[design.block[i]]])
    }
    if(length(design.block) == nblock){
      if(sum(ncomp.design) < ncomp){
        warning("Only design blocks in combination with too high 'ncomp'. Reducing 'ncomp'.")
        ncomp <- sum(ncomp.design)
      }
    }
  }

  # Initialization
  T  <- matrix(0.0, nobj, ncomp) # Scores
  Wb <- list()                   # Block loading weights
  for(i in 1:(nblock+1)){
    Wb[[i]] <- matrix(0.0, m[i], ncomp)
  }
  W <- matrix(0.0, m[nblock+1], ncomp) # Loading weights
  w <- list()
  combo.indexes <- matrix(FALSE, m[nblock+1],ncomb)
  for(i in 1:ncomb){
    u <- unique(combos[i,])
    for(j in 1:length(u)){
      combo.indexes[(1:m[u[j]])+ifelse(u[j]>1, sum(m[1:(u[j]-1)]), 0),i] <- TRUE
    }
  }
  if(!is.null(fixed.order)){ # All blocks fixed
    ncomb <- 1
    fixed <- TRUE
  } else {
    fixed <- FALSE
  }
  if(!is.null(fixed.block)){ # Some blocks fixed
    some.fixed <- TRUE
  } else {
    some.fixed <- FALSE
  }
  tt <- matrix(0.0, nobj, ncomb) # Candidate scores
  tcv <- matrix(0.0, nobj, ncomb) # Candidate scores (cross-validated)
  r  <- array(0.0, dim = c(nobj, nresp, ncomb)) # Projections
  C  <- matrix(0.0, ncomb, ncomp) # Block correlations
  F  <- matrix(0.0, ncomb, ncomp) # Block fit loss
  fitted <- array(0, c(nobj, nresp, ncomp))

  ## Main loop
  for(a in 1:ncomp){
    for(i in 1:ncomb){ # Candidate loading weights and scores
      if(fixed){
        comb <- fixed.order[[a]]
      } else {
        comb <- unique(combos[i,])
      }
      if(some.fixed && length(fixed.block) >= a && !is.null(fixed.block[[a]])){
        comb <- fixed.block[[a]]
      }
      if(canonical){
        if(lower == 0.5 && upper == 0.5){
          # Strategy: Cross-validate t, output also non-cross-validated w and t
          Rlist <- Rcal(X[comb], cbind(y,Y.add), y, weights, nobj, length(comb), internal.validation) 			   # Default CPLS algorithm
        } else {
          Rlist <- RcalP(X[comb], cbind(y,Y.add), y, weights, lower,upper, nobj, length(comb), internal.validation) 			   # Default CPLS algorithm
        }
        w[[i]]  <- Rlist$w
        tt[,i]  <- Rlist$t
        tcv[,i]  <- Rlist$tcv
        # cc[a] <- Rlist$cc
      } else { # PLS(2) type loading weights
        if(length(comb) == 1){
          Xy <- crossprod(X[[comb]], y)
          if(ncol(Xy) > 2 && nrow(Xy) > 2){
            w[[i]]  <- svds(Xy, 1, 1, 0)$u
            tt[,i]  <- do.call(cbind,X[comb]) %*% w[[i]]
            if(!is.logical(internal.validation)){
              nseg <- max(internal.validation)
              for(v in 1:nseg){
                wcv <- svds(crossprod(X[[comb]][internal.validation!=v,], y[internal.validation!=v]))$u
                tcv[internal.validation==v,i] <- X[[comb]][internal.validation!=v,] %*% wcv
              }
            } else {
              tcv[,i] <- tt[,i]
            }
          } else {
            w[[i]]  <- svd(Xy)$u[,1]
            tt[,i]  <- do.call(cbind,X[comb]) %*% w[[i]]
            if(!is.logical(internal.validation)){
              nseg <- max(internal.validation)
              for(v in 1:nseg){
                wcv <- svd(crossprod(X[[comb]][internal.validation!=v,], y[internal.validation!=v]))$u[,1]
                tcv[internal.validation==v,i] <- X[[comb]][internal.validation==v,] %*% wcv
              }
            } else {
              tcv[,i] <- tt[,i]
            }
          }
        } else { # Combined components
          ww <- list(length(comb))
          ttt <- matrix(0, nobj, length(comb))
          for(c in 1:length(comb)){
            ww[[c]] <- svds(crossprod(X[[comb[c]]], y), 1, 1, 0)$u
            ttt[,c] <- X[[comb[c]]] %*% ww[[c]]
          }
          ABr <- cancorr(tt,y,NULL,FALSE)
          w[[i]] <- numeric(0)
          for(c in 1:length(comb)){
            w[[i]] <- c(w[[i]], ABr$A[c,1]*ww[[c]])
          }
          w[[i]] <- w[[i]] / sqrt(c(crossprod(w[[i]])))
          tt[,i]  <- do.call(cbind,X[comb]) %*% w[[i]]
          
          if(!is.logical(internal.validation)){
            nseg <- max(internal.validation)
            for(v in 1:nseg){
              ww <- list(length(comb))
              ttt <- matrix(0, sum(internal.validation!=v), length(comb))
              for(c in 1:length(comb)){
                ww[[c]] <- svds(crossprod(X[[comb[c]]][internal.validation!=v,], y[internal.validation!=v,,drop=FALSE]), 1, 1, 0)$u
                ttt[,c] <- X[[comb[c]]][internal.validation!=v,] %*% ww[[c]]
              }
              ABr <- cancorr(ttt,y,NULL,FALSE)
              wcv <- numeric(0)
              for(c in 1:length(comb)){
                wcv <- c(wcv, ABr$A[c,1]*ww[[c]])
              }
              wcv <- wcv / sqrt(c(crossprod(wcv)))
              tcv[internal.validation==v,i] <- X[[comb[c]]][internal.validation==v,] %*% wcv
            }
          } else {
            tcv[,i] <- tt[,i]
          }
        }
      }

      # Handle design blocks
      if(!is.null(design.block) && any(comb %in% design.block)){
        mdb <- match(comb, design.block); mdb <- mdb[!is.na(mdb)]
        if(any(ncomp.design[mdb] < 1)){
          tt[,i] <- tt[,i] * Inf
          tcv[,i] <- tcv[,i] * Inf
        }
      }
    }

    if(a > 1){ # Orthogonalize scores on previous winning scores
      for(i in 1:ncomb){
        tt[,i] <- tt[,i] - T[, 1:(a-1), drop=FALSE] %*% crossprod(T[, 1:(a-1), drop=FALSE], tt[,i])
        tcv[,i] <- tcv[,i] - T[, 1:(a-1), drop=FALSE] %*% crossprod(T[, 1:(a-1), drop=FALSE], tcv[,i])
      }
    }

    for(i in 1:ncomb){
      tt[,i] <- tt[,i]/sqrt(c(crossprod(tt[,i])))
      tcv[,i] <- tcv[,i]/sqrt(c(crossprod(tcv[,i])))
      yhat   <- tcv[,i] %*% crossprod(tcv[,i],y)
      r[,,i] <- y-yhat
    }

    rsum  <- apply(r^2,3,sum)
    mn    <- min(rsum); mr <- which.min(rsum) # Winner
    if(!fixed){
      C[,a] <- cor(tt[,mr],tt); C[mr,a] <- 1    # Correlation between block scores.
      F[,a] <- sqrt(apply(r^2,3,mean))          # Block-wise fit to residual response.
    }

    # Handle design blocks
    if(!is.null(design.block) && any(unique(combos[mr,] %in% design.block))){
      mdb <- match(unique(combos[mr,]), design.block)
      ncomp.design[mdb] <- ncomp.design[mdb] - 1
    }

    # Handle winners, orthogonalize
    this.fixed <- FALSE
    if(fixed){
      umr <- fixed.order[[a]]
      mr  <- 1
      this.fixed <- TRUE
      combo.fixed <- apply(combo.indexes[,umr,drop=FALSE],1,sum)==1
    } else {
      if(some.fixed && length(fixed.block) >= a && !is.null(fixed.block[[a]])){
        umr <- fixed.block[[a]]
        mri <- NA
        mr  <- 1
        combo.fixed <- apply(combo.indexes[,umr,drop=FALSE],1,sum)==1
        this.fixed <- TRUE
      } else {
        umr <- unique(combos[mr,])
        mri <- mr
      }
    }
    y <- as.matrix(r[,,mr])
    order[[a]] <- umr
    wmr <- rep(0,m[nblock+1])
    if(this.fixed){
      wmr[combo.fixed] <- w[[mr]]
    } else {
      wmr[combo.indexes[,mri]] <- w[[mr]]
    }
    wmr <- wmr - W %*% crossprod(W,wmr)
    count[umr] <- count[umr] + 1;
    W[,a] <- wmr
    T[,a] <- tt[,mr]
  }

  # Loop-less calculations
  P   <- crossprod(X.concat, T) # X-loadings
  PtW <- crossprod(P,W); PtW[lower.tri(PtW)] <- 0 # The W-coordinates of (the projected) P.
  R   <- mrdivide(W,PtW)        # The "SIMPLS weights"
  q   <- crossprod(y.orig,T)    # Regression coeffs (Y-loadings) for the orthogonal scores
  U   <- y.orig %*% q / rep(colSums(q^2), each=nresp)
  if(ncomp > 1)
    for(a in 2:ncomp)
      U[,a] <- U[,a] - T[,1:a] %*% crossprod(T[,1:a], U[,a])

  # Coefficients
  beta <- array(0, dim = c(npred, ncomp, nresp))
  for(i in 1:nresp){
    beta[,,i] <- t(apply(R*rep(q[i,], each=npred),1,cumsum)) # The X-regression coefficients
  }
  beta <- aperm(beta,c(1,3,2))

  # Fitted values
  for(a in 1:ncomp){
    fitted[,,a] <- X.concat %*% beta[,,a]
  }
  fitted <- fitted + rep(Ymeans, each = nobj) # Add mean
  residuals <- - fitted + c(y.orig)
  for(i in 1:(nblock+1)){
    Wb[[i]] <- Wb[[i]][,1:count[i]]
  }

  ## Add dimnames:
  dnX <- dimnames(X.orig)
  dnY <- dimnames(y.orig)
  objnames <- dnX[[1]]
  if (is.null(objnames)) objnames <- dnY[[1]]
  prednames <- dnX[[2]]
  respnames <- dnY[[2]]
  compnames <- paste("Comp", 1:ncomp)
  nCompnames <- paste(1:ncomp, "comps")
  dimnames(T) <- dimnames(U) <- list(objnames, compnames)
  dimnames(W) <- dimnames(P) <-
    list(prednames, compnames)
  dimnames(q) <- list(respnames, compnames)
  dimnames(beta) <- list(prednames, respnames, nCompnames)
  dimnames(fitted) <- dimnames(residuals) <-
    list(objnames, respnames, nCompnames)
  if(ncomb == nblock)
    dimnames(C) <- dimnames(F) <- list(names(X), compnames)
  # colnames(A) <- compnames
  class(T) <- class(U) <- "scores"
  class(P) <- class(W) <- class(q) <- "loadings"
  
  # Explained variance for X
  Xvar <- colSums(P*P); Xtotvar <- sum(X.concat * X.concat)
  attr(T, 'explvar') <- attr(P, 'explvar') <- attr(W, 'explvar') <- Xvar/Xtotvar*100

  list(coefficients=beta, loading.weights=W, loadings=P, scores=T,
       Yloadings = q, Yscores = U, projection=R, PtW=PtW, block.loadings=Wb,
       order=order, count=count, candidate.correlation=C, candidate.RMSE=F,
       Xmeans=Xmeans, Ymeans=Ymeans,
       ncomp=ncomp, X.concat=X.concat, X.orig = X.orig,
       Xvar=Xvar, Xtotvar = Xtotvar,
       fitted.values = fitted, residuals = residuals)
}

#######################
## Rcal function (CPLS)
Rcal <- function(X, Y, Yprim, weights, nobj, nblock, internal.validation) {
  # X <- do.call(cbind,X)
  # W0 <- crossprod(X,Y)
  # Ar <- cancorr(X%*%W0, Yprim, weights, FALSE) # Computes canonical correlations between columns in XW and Y with rows weighted according to 'weights'
  # w  <- W0 %*% Ar$A[,1, drop=FALSE]  # Optimal loadings
  # ifelse(exists('Ar'), a <- Ar$A[,1], a <- NA)
  # list(w = w, cc = Ar$r^2, a = a)

  # # Løsning med Z = [X1*W10 X2*W20 ...]
  Z  <- matrix(0, nobj, nblock*ncol(Y))
  W0 <- list()
  for(b in 1:nblock){
    W0[[b]] <- crossprod(X[[b]],Y)
    Z[ ,ncol(Y)*(b-1)+(1:ncol(Y))] <- X[[b]]%*%W0[[b]]
  }
  Ar <- cancorr(Z, Yprim, weights, FALSE) # Computes canonical correlations between columns in XW and Y with rows weighted according to 'weights'
  w  <- list()
  for(b in 1:nblock){
    w[[b]]  <- W0[[b]] %*% Ar$A[ncol(Y)*(b-1)+(1:ncol(Y)),1, drop=FALSE]  # Optimal loadings
  }
  w <- unlist(w)
  w <- w/sqrt(c(crossprod(w)))
  t <- do.call(cbind, X) %*% w
  tcv <- t
  ifelse(exists('Ar'), a <- Ar$A[,1], a <- NA)
  
  # Cross-validation
  if(!is.logical(internal.validation)){
    nseg <- max(internal.validation)
    for(i in 1:nseg){
      Z  <- matrix(0, sum(internal.validation!=i), nblock*ncol(Y))
      W0 <- list()
      for(b in 1:nblock){
        W0[[b]] <- crossprod(X[[b]][internal.validation!=i,],Y[internal.validation!=i,])
        Z[ ,ncol(Y)*(b-1)+(1:ncol(Y))] <- X[[b]][internal.validation!=i,]%*%W0[[b]]
      }
      Ar <- cancorr(Z, Yprim[internal.validation!=i,], weights[internal.validation!=i], FALSE) # Computes canonical correlations between columns in XW and Y with rows weighted according to 'weights'
      wcv  <- list()
      for(b in 1:nblock){
        wcv[[b]]  <- W0[[b]] %*% Ar$A[ncol(Y)*(b-1)+(1:ncol(Y)),1, drop=FALSE]  # Optimal loadings
      }
      wcv <- unlist(wcv)
      wcv <- wcv/sqrt(c(crossprod(wcv)))
      tcv[internal.validation==i,] <- do.call(cbind, X)[internal.validation==i,,drop=FALSE] %*% wcv
    }
  }
  list(w = w, t = t, tcv = tcv, cc = Ar$r^2, a = a)
}

#########################
## RcalP function (CPPLS)
RcalP <- function(X, Y, Yprim, weights, lower, upper, nobj, nblock, internal.validation){
  lw <- C <- sng <- S <- list()
  for(b in 1:nblock){
    CS <- CorrXY(X[[b]], Y, weights)     # Matrix of corr(Xj,Yg) and vector of std(Xj)
    sng[[b]] <- sign(CS$C)               # Signs of C {-1,0,1}
    C[[b]]  <- abs(CS$C)                 # Correlation without signs
    mS <- max(CS$S); S[[b]] <- CS$S / mS # Divide by largest value
    mC <- max(C[[b]]); C[[b]] <- C[[b]] / mC       #  -------- || --------

    ## Computation of the best vector of loadings
    lw[[b]] <- lw_bestpar(X[[b]], S[[b]], C[[b]], sng[[b]], Yprim, weights, lower, upper, FALSE)
  }
  if(!is.logical(internal.validation)){
    lwcv <- Ccv <- sngcv <- Scv <- list()
    nseg <- max(internal.validation)
    for(i in 1:nseg){
      lwcv[[i]] <- Ccv[[i]] <- sngcv[[i]] <- Scv[[i]] <- list()
      for(b in 1:nblock){
        CS <- CorrXY(X[[b]][internal.validation!=i,], Y[internal.validation!=i,,drop=FALSE], weights[internal.validation!=i])     # Matrix of corr(Xj,Yg) and vector of std(Xj)
        sngcv[[i]][[b]] <- sign(CS$C)               # Signs of C {-1,0,1}
        Ccv[[i]][[b]]  <- abs(CS$C)                 # Correlation without signs
        mS <- max(CS$S); Scv[[i]][[b]] <- CS$S / mS # Divide by largest value
        mC <- max(Ccv[[i]][[b]]); Ccv[[i]][[b]] <- Ccv[[i]][[b]] / mC       #  -------- || --------
        
        ## Computation of the best vector of loadings
        lwcv[[i]][[b]] <- lw_bestpar(X[[b]][internal.validation!=i,], Scv[[i]][[b]], Ccv[[i]][[b]], sngcv[[i]][[b]], Yprim[internal.validation!=i,,drop=FALSE], weights[internal.validation!=i], lower, upper, FALSE)
      }
    }
  }
  if(length(X)==1){ # Single block
    w <- lw[[1]]$w
    w <- w/sqrt(c(crossprod(w)))
    t <- X[[1]] %*% w
    tcv <- t
    if(!is.logical(internal.validation)){
      for(i in 1:nseg){
        wcv <- lwcv[[i]][[1]]$w
        wcv <- wcv/sqrt(c(crossprod(wcv)))
        tcv[internal.validation==i,] <- X[[1]][internal.validation==i,] %*% wcv
      }
    }
    return(list(w=w, t=t, tcv=tcv))

  } else { # Combined blocks
    nresp <- dim(Y)[2]
    W0 <- list()
    Z0 <- matrix(0,nobj,nresp*length(X))
    for(b in 1:nblock){
      if (lw[[b]]$pot == 0) {        # Variable selection from standard deviation
        S[[b]][S[[b]] < max(S[[b]])] <- 0
        W0[[b]] <- S[[b]]
      } else if (lw[[b]]$pot == 1) { # Variable selection from correlation
        C[[b]][C[[b]] < max(C[[b]])] <- 0
        W0[[b]] <- rowSums(C[[b]])
      } else {                     # Standard deviation and correlation with powers
        p <- lw[[b]]$pot                  # Power from optimization
        S[[b]] <- S[[b]]^((1-p)/p)
        W0[[b]] <- (sng[[b]]*(C[[b]]^(p/(1-p))))*S[[b]]
      }
      Z <- X[[b]] %*% W0[[b]]                   # Transform X into W
      Z0[,(1:nresp)+(b-1)*nresp] <- Z
    }
    Ar <- cancorr(Z0, Yprim, weights, FALSE) # Computes canonical correlations between columns in XW and Y with rows weighted according to 'weights'
    w  <- list()
    for(b in 1:nblock){
      w[[b]]  <- W0[[b]] %*% Ar$A[ncol(Y)*(b-1)+(1:ncol(Y)),1, drop=FALSE]  # Optimal loadings
    }
    ifelse(exists('Ar'), a <- Ar$A[,1], a <- NA)
    w <- unlist(w)
    w <- w/sqrt(c(crossprod(w)))
    t <- do.call(cbind, X) %*% w
    tcv <- t
    if(!is.logical(internal.validation)){
      for(i in 1:nseg){
        W0 <- list()
        Z0 <- matrix(0,sum(internal.validation!=i),nresp*length(X))
        for(b in 1:nblock){
          if (lwcv[[i]][[b]]$pot == 0) {        # Variable selection from standard deviation
            Scv[[i]][[b]][S[[b]] < max(Scv[[i]][[b]])] <- 0
            W0[[b]] <- Scv[[i]][[b]]
          } else if (lwcv[[i]][[b]]$pot == 1) { # Variable selection from correlation
            Ccv[[i]][[b]][Ccv[[i]][[b]] < max(Ccv[[i]][[b]])] <- 0
            W0[[b]] <- rowSums(Ccv[[i]][[b]])
          } else {                     # Standard deviation and correlation with powers
            p <- lwcv[[i]][[b]]$pot                  # Power from optimization
            Scv[[i]][[b]] <- Scv[[i]][[b]]^((1-p)/p)
            W0[[b]] <- (sngcv[[i]][[b]]*(Ccv[[i]][[b]]^(p/(1-p))))*Scv[[i]][[b]]
          }
          Z <- X[[b]][internal.validation!=i,] %*% W0[[b]]                   # Transform X into W
          Z0[,(1:nresp)+(b-1)*nresp] <- Z
        }
        Ar <- cancorr(Z0, Yprim[internal.validation!=i,], weights[internal.validation!=i], FALSE) # Computes canonical correlations between columns in XW and Y with rows weighted according to 'weights'
        wcv  <- list()
        for(b in 1:nblock){
          wcv[[b]]  <- W0[[b]] %*% Ar$A[ncol(Y)*(b-1)+(1:ncol(Y)),1, drop=FALSE]  # Optimal loadings
        }
        wcv <- unlist(wcv)
        wcv <- wcv/sqrt(c(crossprod(wcv)))
        tcv[internal.validation==i] <- do.call(cbind, X)[internal.validation==i,] %*% wcv
      }
    }
    return(list(w = w, t = t, tcv = tcv, cc = Ar$r^2, a = a))
  }
}


################
## lw_bestpar function
lw_bestpar <- function(X, S, C, sng, Yprim, weights, lower, upper, trunc.pow) {
  if(!is.null(weights))
    weights <- sqrt(weights) # Prepare weights for cca
  # Compute for S and each columns of C the distance from the median scaled to [0,1]
  if(trunc.pow){
    medC <- t(abs(t(sng*C)-apply(sng*C,2,median)))
    medC <- t(t(medC)/apply(medC,2,max))
    medS <- abs(S-median(S))
    medS <- medS/max(medS)
  } else {
    medS <- medC <- NULL
  }

  #########################
  # Optimization function #
  #########################
  f <- function(p, X, S, C, sng, Yprim, weights, trunc.pow, medS, medC) {
    if(p == 0){         # Variable selection from standard deviation
      S[S < max(S)] <- 0
      W0 <- S
    } else if(p == 1) { # Variable selection from correlation
      C[C < max(C)] <- 0
      W0 <- rowSums(C)
    } else {            # Standard deviation and correlation with powers
      if(trunc.pow){
        ps <- (1-p)/p
        if(ps<1){
          S <- S^ps
        } else {
          S[medS<(1-2*p)] <- 0
        }
        pc <- p/(1-p)
        if(pc<1){
          W0 <- (sng*(C^pc))*S
        } else {
          C[medC<(2*p-1)] <- 0
          W0 <- (sng*C)*S
        }
      } else {
        S <- S^((1-p)/p)
        W0 <- (sng*(C^(p/(1-p))))*S
      }
    }
    Z <- X %*% W0  # Transform X into W0
    -(cancorr(Z, Yprim, weights))^2
  }

  #####################################
  # Logic for optimization segment(s) #
  #####################################
  nOpt <- length(lower)
  pot  <- numeric(3*nOpt)
  ca   <- numeric(3*nOpt)

  for (i in 1:nOpt){
    ca[1+(i-1)*3]  <- f(lower[i], X, S, C, sng, Yprim, weights, trunc.pow, medS, medC)
    pot[1+(i-1)*3] <- lower[i]
    if (lower[i] != upper[i]) {
      Pc <- optimize(f = f, interval = c(lower[i], upper[i]),
                     tol = 10^-4, maximum = FALSE,
                     X = X, S = S, C = C, sng = sng, Yprim = Yprim,
                     weights = weights, trunc.pow = trunc.pow, medS, medC)
      pot[2+(i-1)*3] <- Pc[[1]]; ca[2+(i-1)*3] <- Pc[[2]]
    }
    ca[3+(i-1)*3]  <- f(upper[i], X, S, C, sng, Yprim, weights, trunc.pow, medS, medC)
    pot[3+(i-1)*3] <- upper[i]
  }


  ########################################################
  # Computation of final w-vectors based on optimization #
  ########################################################
  cc <- max(-ca)                      # Determine which is more succesful
  cmin <- which.max(-ca)              # Determine which is more succesful
  if (pot[cmin] == 0) {        # Variable selection from standard deviation
    S[S < max(S)] <- 0
    w <- S
  } else if (pot[cmin] == 1) { # Variable selection from correlation
    C[C < max(C)] <- 0
    w <- rowSums(C)
  } else {                     # Standard deviation and correlation with powers
    p <- pot[cmin]                  # Power from optimization
    if(trunc.pow){ # New power algorithm
      ps <- (1-p)/p
      if(ps<1){
        S <- S^ps
      } else {
        S[medS<(1-2*p)] <- 0
      }
      pc <- p/(1-p)
      if(pc<1){
        W0 <- (sng*(C^pc))*S
      } else {
        C[medC<(2*p-1)] <- 0
        W0 <- (sng*C)*S
      }
    } else {
      S <- S^((1-p)/p)
      W0 <- (sng*(C^(p/(1-p))))*S
    }

    Z <- X %*% W0                   # Transform X into W
    Ar <- cancorr(Z, Yprim, weights, FALSE) # Computes canonical correlations between columns in XW and Y with rows weighted according to 'weights'
    w <- W0 %*% Ar$A[,1, drop=FALSE]  # Optimal loadings
  }
  pot <- pot[cmin]
  ifelse(exists('Ar'), a <- Ar$A[,1], a <- NA)
  list(w = w, pot = pot, cc = cc, a = a)
}


################
## CorrXY function
CorrXY <- function(X, Y, weights) {
  ##  Computation of correlations between the columns of X and Y
  n  <- dim(X)[1]
  if(is.null(weights)){
    cx <- colMeans(X)
    cy <- colMeans(Y)
    X <- X - rep(cx, each = n)
    Y <- Y - rep(cy, each = n)
  } else {
    cx <- crossprod(weights,X)/sum(weights)
    cy <- crossprod(weights,Y)/sum(weights)
    X  <- X - rep(cx, each = n)
    Y  <- Y - rep(cy, each = n)
    X  <- X * weights
    Y  <- Y * weights
  }

  sdX <- sqrt(apply(X^2,2,mean))
  inds <- which(sdX == 0, arr.ind=FALSE)
  sdX[inds] <- 1

  ccxy <- crossprod(X, Y) / (n * tcrossprod(sdX, sqrt(apply(Y^2,2,mean))))
  sdX[inds] <- 0
  ccxy[inds,] <- 0
  list(C = ccxy, S = sdX)
}

################
## function norm
#norm <- function(vec) {
#  sqrt(crossprod(vec)[1])
#}

################
## Stripped version of canonical correlation (cancor)
cancorr <- function (x, y, weights, opt = TRUE) {
  nr  <- nrow(x)
  ncx <- ncol(x)
  ncy <- ncol(y)
  if (!is.null(weights)){
    x <- x * weights
    y <- y * weights
  }
  qx <- qr(x, LAPACK = TRUE)
  qy <- qr(y, LAPACK = TRUE)
  qxR <- qr.R(qx)
  # Compute rank like MATLAB does
  dx <- sum(abs(diag(qxR)) > .Machine$double.eps*2^floor(log2(abs(qxR[1])))*max(nr,ncx))
  if (!dx)
    stop("'x' has rank 0")
  qyR <- qr.R(qy)
  # Compute rank like MATLAB does
  dy <- sum(abs(diag(qyR)) > .Machine$double.eps*2^floor(log2(abs(qyR[1])))*max(nr,ncy))
  if (!dy)
    stop("'y' has rank 0")
  dxy <- min(dx,dy)
  if(opt) {
    z <- svd(qr.qty(qx, qr.qy(qy, diag(1, nr, dy)))[1:dx,, drop = FALSE],
             nu = 0, nv = 0)
    ret <- max(min(z$d[1],1),0)
  } else {
    z <- svd(qr.qty(qx, qr.qy(qy, diag(1, nr, dy)))[1:dx,, drop = FALSE],
             nu = dxy, nv = 0)
    A <- backsolve((qx$qr)[1:dx,1:dx, drop = FALSE], z$u)*sqrt(nr-1)
    if((ncx - nrow(A)) > 0) {
      A <- rbind(A, matrix(0, ncx - nrow(A), dxy))
    }
    A[qx$pivot,] <- A
    ret <- list(A=A, r=max(min(z$d[1],1),0))
  }
  ret
}

# All possible sorted block combinations up to nbmax length
common.comb <- function(nb, nbmax){
  if(nbmax == 1)
    return(matrix(1:nb, ncol=1))
  comb <- list()
  for(i in 1:nbmax){
    comb[[i]] <- 1:nb
  }
  comb  <- unique(t(apply(do.call(expand.grid,comb)[,nbmax:1],1,sort)))
  diffs <- apply(comb,1,diff)
  rbind(comb[diffs==0,],comb[diffs!=0,])
}
