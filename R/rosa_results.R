#' @name rosa_results
#' @title Result functions for ROSA models
#'
#' @aliases predict.rosa rosa.classify coef.rosa print.rosa summary.rosa blockexpl print.rosaexpl scores.rosa loadings.rosa
#' @param object A \code{rosa} object.
#' @param x A \code{rosa} object.
#' @param newdata Optional new data with the same types of predictor blocks as the ones used for fitting the object.
#' @param ncomp An \code{integer} giving the number of components to apply (cummulative).
#' @param comps An \code{integer} vector giving the exact components to apply (subset).
#' @param type For \code{predict}: A \code{character} indicating if responses or scores should be predicted (default = "response", or "scores")
#' @param type For \code{blockexpl}: Character indicating which type of explained variance to compute (default = "train", alternative = "CV").
#' @param na.action Function determining what to do with missing values in \code{newdata}.
#' @param intercept A \code{logical} indicating if coefficients for the intercept should be included (default = FALSE).
#' @param what A \code{character} indicating if summary should include all, validation or training.
#' @param digits The number of digits used for printing.
#' @param print.gap Gap between columns when printing.
#' @param classes A \code{character} vector of class labels.
#' @param LQ A \code{character} indicating if 'max' (maximum score value), 'lda' or 'qda' should be used when classifying.
#' @param compwise Logical indicating if block-wise (default/FALSE) or component-wise (TRUE) explained variance should be printed.
#' @param ... Additional arguments. Currently not implemented.
#'
#' @return Returns depend on method used, e.g. \code{predict.rosa} returns predicted responses 
#' or scores depending on inputs, \code{coef.rosa} returns regression coefficients, \code{blockexpl}
#' returns an object of class \code{rosaexpl} containing block-wise and component-wise explained variance contained in a matrix with attributes.
#' 
#' @description Standard result computation and extraction functions for ROSA (\code{\link{rosa}}).
#' 
#' @details Usage of the functions are shown using generics in the examples below.
#' Prediction, regression coefficients, object printing and summary are available through: 
#' \code{predict.rosa},  \code{coef.rosa}, \code{print.rosa} and \code{summary.rosa}.
#' Explained variances are available (block-wise and global) through \code{blockexpl} and \code{print.rosaexpl}.
#' Scores and loadings have their own extensions of \code{scores()} and \code{loadings()} throught
#' \code{scores.rosa} and \code{loadings.rosa}. Finally, there is work in progress on classifcation
#' support through \code{rosa.classify}.
#' 
#' When \code{type} is \code{"response"} (default), predicted response values
#' are returned.  If \code{comps} is missing (or is \code{NULL}), predictions
#' for \code{length(ncomp)} models with \code{ncomp[1]} components,
#' \code{ncomp[2]} components, etc., are returned.  Otherwise, predictions for
#' a single model with the exact components in \code{comps} are returned.
#' (Note that in both cases, the intercept is always included in the
#' predictions.  It can be removed by subtracting the \code{Ymeans} component
#' of the fitted model.)
#' 
#' If \code{comps} is missing (or is \code{NULL}), \code{coef()[,,ncomp[i]]}
#' are the coefficients for models with \code{ncomp[i]} components, for \eqn{i
#' = 1, \ldots, length(ncomp)}.  Also, if \code{intercept = TRUE}, the first
#' dimension is \eqn{nxvar + 1}, with the intercept coefficients as the first
#' row.
#'
#' If \code{comps} is given, however, \code{coef()[,,comps[i]]} are the
#' coefficients for a model with only the component \code{comps[i]}, i.e., the
#' contribution of the component \code{comps[i]} on the regression
#' coefficients.
#' 
#' @references Liland, K.H., Næs, T., and Indahl, U.G. (2016). ROSA - a fast extension of partial least squares regression for multiblock data analysis. Journal of Chemometrics, 30, 651–662, doi:10.1002/cem.2824.
#'
#' @examples
#' data(potato)
#' mod <- rosa(Sensory[,1] ~ ., data = potato, ncomp = 5, subset = 1:20)
#' testset <- potato[-(1:20),]; testset$Sensory <- testset$Sensory[,1,drop=FALSE]
#' predict(mod, testset, ncomp=5)
#' dim(coef(mod, ncomp=5)) # <variables x responses x components>
#' print(mod)
#' summary(mod)
#' blockexpl(mod)
#' print(blockexpl(mod), compwise=TRUE)
#' 
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results and plotting are found in \code{\link{rosa_results}} and \code{\link{rosa_plots}}, respectively.
#' @export
predict.rosa <- function(object, newdata, ncomp = 1:object$ncomp, comps,
                         type = c("response", "scores"), na.action = na.pass, ...){
  if (missing(newdata) || is.null(newdata)){
    # newX <- object$X.concat
    newX <- model.matrix(object)
    colnames(newX) <- colnames(object$X.concat)
  } else {
    if(is.matrix(newdata)){
      if (ncol(newdata) != length(object$Xmeans))
        stop("'newdata' does not have the correct number of columns")
      newX <- newdata
    } else { # Assume newdata is a list
      newX <- model.frame(formula(object), data = newdata)
      newX <- do.call(cbind,newX[-1])
      # Check for missing dimnames
      if(is.null(rownames(newX)))
        rownames(newX) <- 1:nrow(newX)
      if(is.null(colnames(newX)))
        colnames(newX) <- 1:ncol(newX)
      if(ncol(newX) != length(object$Xmeans))
        stop("'newdata' does not have the correct number of columns")
    }
  }

  nobs <- dim(newX)[1]

  if (!is.null(object$scale)) newX <- newX / rep(object$scale, each = nobs)
  type <- match.arg(type)

  if (type == "response") {
    if (missing(comps) || is.null(comps)) {
      ## Predict with models containing ncomp[1] components,
      ## ncomp[2] components, etc.
      if (missing(newdata)) return(fitted(object)[,,ncomp, drop=FALSE])
      B <- coef(object, ncomp = ncomp, intercept = TRUE)
      dPred <- dim(B)
      dPred[1] <- dim(newX)[1]
      dnPred <- dimnames(B)
      dnPred[1] <-
        if(is.null(dimnames(newX))) list(NULL) else dimnames(newX)[1]
      pred <- array(dim = dPred, dimnames = dnPred)
      for (i in seq(along = ncomp))
        pred[,,i] <- newX %*% B[-1,,i] + rep(B[1,,i], each = nobs)
      return(pred)
    } else {
      ## Predict with a model containing the components `comps'
      B <- rowSums(coef(object, comps = comps), dims = 2)
      B0 <- object$Ymeans - object$Xmeans %*% B
      pred <- newX %*% B + rep(c(B0), each = nobs)
      if (missing(newdata) && !is.null(object$na.action))
        pred <- napredict(object$na.action, pred)
      return(pred)
    }
  } else {
    ## Return predicted scores (for scores, `cumulative' has no meaning)
    ## When predicting scores, we allow ncomp as an alias for comps:
    if (missing(comps) || is.null(comps)) comps <- ncomp
    if (missing(newdata)) {
      TT <- object$scores[,comps]
      if (!is.null(object$na.action))  TT <- napredict(object$na.action, TT)
    } else {
      if (is.null(object$projection))
        stop("`object' has no `projection' component.  Maybe it was fitted with `stripped = TRUE'.")
      TT <- (newX - rep(object$Xmeans, each = nobs)) %*%
        object$projection[,comps]
    }
    return(TT)
  }
}

#' @rdname rosa_results
#' @export
coef.rosa <- function(object, ncomp = object$ncomp, comps, intercept = FALSE,
                     ...)
{
  if (missing(comps) || is.null(comps)) {
    ## Cumulative coefficients:
    B <- object$coefficients[,,ncomp, drop=FALSE]
    if (intercept == TRUE) {      # Intercept has only meaning for
      # cumulative coefficients
      dB <- dim(B)
      dB[1] <- dB[1] + 1
      dnB <- dimnames(B)
      dnB[[1]] <- c("(Intercept)", dnB[[1]])
      BInt <- array(dim = dB, dimnames = dnB)
      BInt[-1,,] <- B
      for (i in seq(along = ncomp))
        BInt[1,,i] <- object$Ymeans - object$Xmeans %*% B[,,i]
      B <- BInt
    }
  } else {
    ## Individual coefficients:
    B <- object$coefficients[,,comps, drop=FALSE]
    g1 <- which(comps > 1)
    ## Indiv. coef. must be calculated since object$coefficients is
    ## cumulative coefs.
    B[,,g1] <- B[,,g1, drop=FALSE] -
      object$coefficients[,,comps[g1] - 1, drop=FALSE]
    dimnames(B)[[3]] <- paste("Comp", comps)
  }
  return(B)
}

#' @rdname rosa_results
#' @export
print.rosa <- function(x, ...) {
     ana <- "Response Orinented Sequential Alternation"
     alg <- "CPPLS"
  cat(ana, ", fitted with the", alg, "algorithm.")
  if (!is.null(x$validation))
    cat("\nCross-validated using", length(x$validation$segments),
        attr(x$validation$segments, "type"), "segments.")
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  invisible(x)
}

#' @rdname rosa_results
#' @export
summary.rosa <- function(object, what = c("all", "validation", "training"),
                        digits = 4, print.gap = 2, ...)
{
  what <- match.arg(what)
  if (what == "all") what <- c("validation", "training")
  if (is.null(object$validation)) what <- "training"

  nobj <- nrow(object$scores)
  nresp <- length(object$Ymeans)
  yvarnames <- dimnames(fitted(object))[[2]]
  cat("Data: \tX dimension:", nobj, length(object$Xmeans),
      "\n\tY dimension:", nobj, nresp)
  cat("\nFit method:", object$method)
  cat("\nNumber of components considered:", object$ncomp)

  for (wh in what) {
    if (wh == "training") {
      cat("\nTRAINING: % variance explained\n")
      xve <- explvar(object)
      yve <- 100 * drop(R2(object, estimate = "train",
                           intercept = FALSE)$val)
      tbl <- rbind(cumsum(xve), yve)
      dimnames(tbl) <- list(c("X", yvarnames),
                            paste(1:object$ncomp, "comps"))
      print(tbl, digits = digits, print.gap = print.gap, ...)
    } else {
      cat("\n\nVALIDATION: RMSEP")
      cat("\nCross-validated using", length(object$validation$segments),
          attr(object$validation$segments, "type"), "segments.\n")
      print(RMSEP(object), digits = digits, print.gap = print.gap, ...)
    }
  }
}

#' @rdname rosa_results
#' @export
blockexpl <- function(object, ncomp = object$ncomp, type = c("train","CV")){
  nblock  <- length(object$X)
  X       <- object$X.orig
  categ   <- !is.null(object$classes)
  Y <- model.response(object$model)
  nobj    <- ifelse(is.null(dim(Y)), length(Y),nrow(Y))
  if(categ){
    Y <- model.matrix(~Y-1, data.frame(Y = factor(Y)))
  }
  Ytotvar <- sum(c(Y-rep(object$Ymeans, each = nobj))^2)
  Xtotvar <- object$Xtotvar
  Xvar    <- object$Xvar

  # Component-wise explained variance
  if(type[1] == "train"){
    fits <- object$fitted.values
  } else {
    fits <- object$validation$pred
    ncomp <- min(ncomp, dim(fits)[3])
  }
  resc <- matrix(0.0, 2, ncomp+1)
  resc[1,1] <- Xvar[1]/Xtotvar
  resc[2,1] <- 1-sum(c(fits[,,1]-Y)^2)/Ytotvar
  for(i in 2:ncomp){
    resc[1,i] <- Xvar[i]/Xtotvar
    resc[2,i] <- 1-sum(c(fits[,,i]-Y)^2)/Ytotvar - sum(resc[2,])
    # resc[2,i] <- sum(c(fits[,,i]-fits[,,i-1])^2)/Ytotvar
  }
  resc[,ncomp+1] <- sum(c(fits[,,i]-Y)^2)/Ytotvar #1-rowSums(resc) # res[,nblock+1]
  nameVec <- character(ncomp)
  for(i in 1:ncomp){
    nameVec[i] <- paste(names(object$X)[object$order[[i]]], collapse = ":")
  }
  dimnames(resc) <- list(c("X","Y"), c(paste(paste("comp.", 1:ncomp, sep=""), paste(' (',nameVec,')',sep=""),sep=""),"residual"))

  cblocks <- unlist(lapply(object$order, function(i)length(i)>1)) # Common blocks
  cnblock <- length(un <- unique(object$order[cblocks]))+nblock
  corder  <- as.list(c(1:nblock, un))

  # Block-wise explained variance
  res     <- matrix(0.0, 2, cnblock+1)
  for(i in 1:cnblock){
    ids <- unlist(lapply(object$order[1:ncomp], function(j)identical(j,corder[[i]])))
    if(sum(ids) > 0){
      res[1,i] <- sum(Xvar[ids]) / Xtotvar #- res[1,1]
      res[2,i] <- sum(resc[,-(ncomp+1)][2,ids])
    }
  }
  if(any(res[!is.na(res)] < 0)){
    warning('Negative block-wise explained variance encountered.')
  }
  res[,cnblock+1] <- 1-rowSums(res, na.rm = TRUE)
  nameVec <- character(cnblock)
  for(i in 1:cnblock){
    nameVec[i] <- paste(names(object$X)[corder[[i]]], collapse=":")
  }
  dimnames(res) <- list(c("X","Y"), c(nameVec,"residual"))
  if(length(remove <- nblock + which(colSums(res[,-c(1:nblock,cnblock+1), drop=FALSE])==0)) > 0){
    res <- res[,-remove]
    corder <- corder[-remove]
  }

  if(type[1] != "train"){
    res  <- res[2,,drop=FALSE]
    resc <- resc[2,,drop=FALSE]
  }

  attr(res,'compwise') <- resc
  attr(res,'index') <- corder
  class(res) <- c("rosaexpl","matrix")
  res
}

#' @rdname rosa_results
#' @export
print.rosaexpl <- function(x, digits = 3, compwise = FALSE, ...){
  if(compwise){
    cat("Component-wise explained variance\n\n")
    print(round(attr(x, "compwise"), digits))
  } else {
    attr(x, "compwise") <- NULL
    attr(x, "index") <- NULL
    class(x) <- "matrix"
    cat("Block-wise explained variance\n\n")
    print(round(x, digits))
  }
}

# TODO: S3 object of rosa.classify #' @method pls ER
#' @rdname rosa_results
#' @importFrom MASS lda qda
#' @export
rosa.classify <- function(object, classes, newdata, ncomp, LQ){
  if(LQ == "max"){
    labels  <- names(table(classes))
    predVal <- predict(object, newdata = newdata, ncomp = 1:ncomp)
    class   <- apply(predVal,c(1,3),which.max)
    for(i in 1:ncol(class)){
      class[[i]]   <- labels[class[[i]]]
    }
    colnames(class) <- paste("Comp.", 1:ncomp, sep="")
    return(class)
    
  } else { # LDA or QDA
    # Extract and predict scores
    scoresCal <- scores(object)
    scoresVal <- predict(object, newdata = newdata, type = "scores")
    
    # Prepare for storage
    N <- dim(scoresVal)
    class <- matrix(0, N[1],ncomp)
    
    # Create ncomp lda models and predict classes
    for(i in 1:ncomp){
      if(LQ == "lda"){
        ldai <- lda(scoresCal[, 1:i, drop = FALSE], classes, tol = 1.0e-10)
      }
      if(LQ == "qda"){
        ldai <- qda(scoresCal[, 1:i, drop = FALSE], classes, tol = 1.0e-10)
      }
      class[, i] <- predict(ldai, scoresVal[, 1:i, drop = FALSE])$class
    }
    colnames(class) <- paste("Comp.", 1:ncomp, sep="")
    return(class)
  }
}

#' @rdname rosa_results
#' @export
scores.rosa <- function(object, ...) {
  object$scores
}

#' @rdname rosa_results
#' @export
loadings.rosa <- function(object, ...) {
  object$loadings
}


################
# Undocumented #
################

# model.matrix.rosa <- function (object, ...) {
#   if (n_match <- match("X.concat", names(object), 0))
#     object[[n_match]]
#   else {
#     data <- model.frame(object, ...)
#     mm <- NextMethod("model.matrix", data = data)
#     mm <- delete.intercept(mm)
#     mt <- terms(object)
#     if (length(attr(mt, "term.labels")) == 1 && !is.null(colnames(data[[attr(mt,
#                                                                              "term.labels")]])))
#       colnames(mm) <- sub(attr(mt, "term.labels"), "",
#                           colnames(mm))
#     return(mm)
#   }
# }


#
# ## Print method for mvrVal objects:
# print.mvrVal <- function(x, digits = 4, print.gap = 2, ...) {
#   nresp <- dim(x$val)[2]
#   yvarnames <- dimnames(x$val)[[2]]
#   names(dimnames(x$val)) <- NULL
#   for (i in 1:nresp) {
#     if (nresp > 1) cat("\nResponse:", yvarnames[i], "\n")
#     print(x$val[,i,], digits = digits, print.gap = print.gap, ...)
#   }
#   invisible(x)
# }
