rosaCV <- function(X, Y, Y.cat, ncomp, response.type, Y.add = NULL, common.comp,
                  scale = FALSE, weights = weights, segments = 10,
                  segment.type = c("random", "consecutive", "interleaved"), fixed.order, canonical,
                  length.seg, jackknife = FALSE, trace = FALSE, ...)
{
  ## Initialise:
  Y <- as.matrix(Y)
  if(!(missing(Y.add) || is.null(Y.add)))
    Y.add <- as.matrix(Y.add)
  ## Save dimnames:
  dnX <- rownames(X[[1]])
  dnY <- dimnames(Y)
  ## Remove dimnames for performance (doesn't seem to matter; in fact,
  ## as far as it has any effect, it hurts a tiny bit in most situations).
  ## dimnames(X) <- dimnames(Y) <- NULL

  nobj <- dim(X[[1]])[1]
  npred <- dim(X[[1]])[2]
  nresp <- dim(Y)[2]
  nblock <- length(X)

  ## Check the `scale' parameter:
  if (!is.logical(scale) || length(scale) != 1)
    stop("'scale' must be 'TRUE' or 'FALSE'")

  ## Set up segments:
  if (is.list(segments)) {
    if (is.null(attr(segments, "type")))
      attr(segments, "type") <- "user supplied"
  } else {
    if (missing(length.seg)) {
      segments <- cvsegments(nobj, k = segments, type = segment.type)
    } else {
      segments <- cvsegments(nobj, length.seg = length.seg,
                             type = segment.type)
    }
  }

  ## Reduce ncomp, if neccessary:
  ncomp <- min(ncomp, nobj - max(sapply(segments, length)) - 1)

  ## Variables to save CV results in:
  adj <- matrix(0, nrow = nresp, ncol = ncomp)
  cvPred <- pred <- array(0, dim = c(nobj, nresp, ncomp))
  cvClass <- matrix(0, nobj,ncomp)
  if (jackknife)
    cvCoef <- array(dim = c(npred, nresp, ncomp, length(segments)))
  gammas <- list()

  if (trace) cat("Segment: ")
  for (n.seg in 1:length(segments)) {
    if (trace) cat(n.seg, "")

    ## Set up train data:
    seg <- segments[[n.seg]]
    Xtrain <- list()
    for(b in 1:nblock){
      Xtrain[[b]] <- X[[b]][-seg,,drop=FALSE]
    }
    sdtrain <- list()
    if (scale) {
      ntrain <- nrow(Xtrain[[1]])
      ## This is faster than sd(X), but cannot handle missing values:
      for(b in 1:nblock){
        sdtrain[[b]] <-
          sqrt(colSums((Xtrain[[b]] - rep(colMeans(Xtrain[[b]]), each = ntrain))^2) /
                 (ntrain - 1))
        if (any(abs(sdtrain[[b]]) < .Machine$double.eps^0.5))
          warning("Scaling with (near) zero standard deviation")
        Xtrain[[b]] <- Xtrain[[b]] / rep(sdtrain[[b]], each = ntrain)
      }
    }
    X.concat <- do.call(cbind, Xtrain)

    ## Fit the model:
    fit <- rosa.fit(Xtrain, X.concat, Y[-seg,,drop=FALSE], Y.add[-seg,,drop=FALSE], ncomp, common.comp, weights[-seg], fixed.order, NULL, NULL, canonical, FALSE)

    ## Optionally collect coefficients:
    if (jackknife) cvCoef[,,,n.seg] <- fit$coefficients

    ## Collect gamma-values from CPPLS
    gammas[[n.seg]] <- fit$gammas

    ## Set up test data:
    Xtests <- list()
    for(b in 1:nblock){
      if(scale){
        Xtests[[b]] <- X[[b]] / rep(sdtrain[[b]], each = nobj)
      } else {
        Xtests[[b]] <- X[[b]]
      }
    }
    Xtest <- do.call(cbind, Xtests)
    Xtest <- Xtest - rep(fit$Xmeans, each = nobj)

    ## Predict test data:
    Ymeansrep <- rep(fit$Ymeans, each = nobj)
    for (a in 1:ncomp)
      pred[,,a] <- Xtest %*% fit$coefficients[,,a] + Ymeansrep

    ## Save the cross-validated predictions:
    cvPred[seg,,] <- pred[seg,,, drop=FALSE]
    adj <- adj + length(seg) * colSums((pred - c(Y))^2)

    ## Classification
    if(response.type == "categorical"){
      class(fit) <- c('rosa','mvr')
      Xt <- list()
      for(b in 1:nblock){
        Xt[[b]] <- Xtests[[b]][seg,,drop=FALSE]
      }
      cvClass[seg,] <- rosa.classify(fit, Y.cat[-seg,,drop=FALSE], Xt, ncomp, 'lda')
    }
  }
  if (trace) cat("\n")

  ## Calculate validation statistics:
  PRESS0 <- apply(Y, 2, var) * nobj^2 / (nobj - 1) # FIXME: Only correct for loocv!
  PRESS <- colSums((cvPred - c(Y))^2)

  ## Add dimnames:
  objnames <- dnX
  if (is.null(objnames)) objnames <- dnY[[1]]
  respnames <- dnY[[2]]
  nCompnames <- paste(1:ncomp, "comps")
  names(PRESS0) <- respnames
  dimnames(adj) <- dimnames(PRESS) <-
    list(respnames, nCompnames)
  dimnames(cvPred) <- list(objnames, respnames, nCompnames)
  if (jackknife)
    dimnames(cvCoef) <- list(NULL, respnames, nCompnames,
                             paste("Seg", seq.int(along.with = segments)))

  list(method = "CV", pred = cvPred, coefficients = if (jackknife) cvCoef,
       gammas = gammas, classes = if(response.type == "categorical") cvClass,
       PRESS0 = PRESS0, PRESS = PRESS, adj = adj / nobj^2,
       segments = segments, ncomp = ncomp)
}
