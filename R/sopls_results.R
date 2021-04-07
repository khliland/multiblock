#' @name sopls_object
#' @title Result functions for SO-PLS models
#'
#' @aliases predict.sopls sopls.classify coef.sopls print.sopls summary.sopls
#' @param object A \code{sopls} object.
#' @param x A \code{sopls} object.
#' @param newdata Optional new data with the same types of predictor blocks as the ones used for fitting the object.
#' @param ncomp An \code{integer} giving the number of components to apply.
#' @param comps An \code{integer} vector giving the exact components to apply.
#' @param type A \code{character} for \code{predict} indicating if responses or scores should be predicted (default = "response", or "scores"), for \code{summary} indicating which type of explained variance to compute (default = "train", alternative = "CV").
#' @param na.action Function determining what to do with missing values in \code{newdata}.
#' @param intercept A \code{logical} indicating if coefficients for the intercept should be included (default = FALSE).
#' @param what A \code{character} indicating if summary should include all, validation or training.
#' @param digits The number of digits used for printing.
#' @param print.gap Gap between columns when printing.
#' @param classes A \code{character} vector of class labels.
#' @param LQ A \code{character} indicating if 'max' (maximum score value), 'lda' or 'qda' should be used when classifying.
#' @param block An \code{integer} indicating which block to use.
#' @param estimate A \code{character} indicating if 'train', 'CV' or 'test' results should be displayed.
#' @param individual A \code{logical} indicating if results for individual responses should be displayed.
#' @param ... Additional arguments. Currently not implemented.
#'
#' @return Returns depend on method used, e.g. \code{predict.sopls} returns predicted responses 
#' or scores depending on inputs, \code{coef.sopls} return regression coefficients.
#' 
#' @description Standard result functions for SO-PLS (\code{\link{sopls}}).
#' 
#' @references Jørgensen K, Mevik BH, Næs T. Combining designed experiments with several blocks of spectroscopic data. Chemometr Intell Lab Syst. 2007;88(2): 154–166.
#'
#' @examples
#' data(potato)
#' mod <- sopls(Sensory[,1] ~ ., data = potato[c(1:3,9)], ncomp = 5, subset = 1:20)
#' testset <- potato[-(1:20),]; testset$Sensory <- testset$Sensory[,1,drop=FALSE]
#' predict(mod, testset, comps=c(2,1,2))
#' #dim(coef(mod, ncomp=5)) # <variables x responses x components>
#' print(mod)
#' summary(mod)
#' 
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' @export
predict.sopls <- function(object, newdata, ncomp = object$ncomp, comps = object$ncomp,
                          type = c("response", "scores"), na.action = na.pass, ...){
  if (missing(newdata) || is.null(newdata)){
    newdata <- object$data$X
  } else {
    newdata <- model.frame(formula(object), data = newdata)
    newdata <- newdata[-1]
  }
  
  type <- match.arg(type)
  
  if(!is.character(comps) && sum(comps) > object$max_comp)
    stop(paste0("Selected components (",paste(comps,collapse=","),") is outside range of fitted model."))
  if (type == "response") {
    return(sopls_prediction(object, newdata, comps, FALSE))
  } else {
    return(sopls_prediction(object, newdata, comps, TRUE))
  }
}

#' @rdname sopls_object
#' @export
coef.sopls <- function(object, ncomp = object$ncomp, comps = object$ncomp, intercept = FALSE,
                       ...)
{
  if(sum(comps) > object$max_comp)
    stop(paste0("Selected components (",paste(comps,collapse=","),") is outside range of fitted model."))
  
  X <- object$data$X
  Y <- object$data$Y
  n    <- dim(X[[1]])[1]
  
  nblock  <- length(X)
  nresp   <- dim(Y)[2]
  selComp <- pathComp(comps, object$decomp$compList)
  tot_comp <- nrow(selComp$path)
  
  Cr <- Crval <- 0
  for(i in 1:nblock){
    if(comps[i]>0){
      X[[i]]    <- as.matrix(X[[i]])
      X[[i]]    <- X[[i]]    - rep(colMeans(X[[i]]), each = n)
      Cr    <- Cr    + tcrossprodQ(X[[i]], X[[i]])    %*% object$decomp$Ry[, selComp$hits, drop=FALSE]
    }
  }
  Xconcat <- do.call(cbind, X)
  # XW(P'W)^-1, ie. WB without Q
  no_Q <- crossprod(Xconcat, object$decomp$Ry[,selComp$hits,drop=FALSE]) %*% solve(crossprodQ(object$decomp$T[,selComp$hits,drop=FALSE], Cr))
  
  # Coefficients per response
  B <- array(0, c(ncol(Xconcat), nresp, tot_comp))  
  for(r in 1:nresp){
    B[,r,] <- t(apply(no_Q * object$decomp$Q[r,selComp$hits], 1, cumsum))
  }
  dimnames(B) <- list(colnames(Xconcat), colnames(Y), apply(selComp$path,1,paste,collapse=","))
  return(B)
}

#' @rdname sopls_object
#' @export
print.sopls <- function(x, ...) {
  ana <- "Sequential and Orthogonalized Partial Least Squares"
  alg <- "PKPLS"
  cat(ana, ", fitted with the ", alg, " algorithm.", sep="")
  if (!is.null(x$validation))
    cat("\nCross-validated using", length(x$validation$segments),
        attr(x$validation$segments, "type"), "segments.")
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  invisible(x)
}

#' @rdname sopls_object
#' @export
summary.sopls <- function(object, what = c("all", "validation", "training"),
                          digits = 4, print.gap = 2, ...)
{
  what <- match.arg(what)
  if (what == "all") what <- c("validation", "training")
  if (is.null(object$validation)) what <- "training"
  
  nobj <- nrow(object$decomp$T)
  nresp <- nrow(object$decomp$Q)
  yvarnames <- rownames(object$decomp$Q)
  cat("Data: \tX dimension:", nobj, sum(unlist(lapply(object$decom$X, ncol))),
      "\n\tY dimension:", nobj, nresp)
  cat("\nFit method:", object$method)
  cat("\nNumber of components considered:", object$max_comps)
  
  for (wh in what) {
    if (wh == "training") {
      cat("\nTRAINING: % variance explained\n")
      Xcat <- scale(do.call(cbind, object$data$X), scale=FALSE)
      P <- crossprod(object$decomp$T, Xcat)
      ssxx <- 100 * apply(P^2,1,sum)/sum(Xcat^2)
      xve <- ssxx
      for(i in 1:length(xve)){
        xve[i] <- sum(ssxx[pathComp(object$decomp$compList[i,], object$decomp$compList)$hits])
      }
      ssyy <- 100 * object$decomp$Q^2/apply(as.matrix(object$model[[1]]),2,function(x){x <- x-mean(x); sum(x^2)})
      yve <- ssyy
      for(i in 1:ncol(yve)){
        yve[,i] <- rowSums(ssyy[, pathComp(object$decomp$compList[i,], object$decomp$compList)$hits, drop=FALSE])
      }
      tbl <- rbind(xve, yve)
      rownames(tbl)[1] <- "X"
      print(tbl, digits = digits, print.gap = print.gap, ...)
    } else {
      cat("\n\nVALIDATION: RMSEP")
      cat("\nCross-validated using", length(object$validation$segments),
          attr(object$validation$segments, "type"), "segments.\n")
      print(object$validation$RMSECV, digits = digits, print.gap = print.gap, ...)
    }
  }
}


#' @rdname sopls_object
#' @export
classify <- function(object, ...) UseMethod("classify")

#' @rdname sopls_object
#' @importFrom MASS lda qda
#' @export
classify.sopls <- function(object, classes, newdata, ncomp, LQ = "LDA", ...){
  if(LQ == "max"){
    labels  <- names(table(classes))
    predVal <- predict(object, newdata = newdata, ncomp = 1:ncomp) # Njei. Ikke sånn
    class   <- apply(predVal,c(1,3),which.max)
    for(i in 1:ncol(class)){
      class[[i]]   <- labels[class[[i]]]
    }
    colnames(class) <- paste("Comp.", 1:ncomp, sep="")
    return(class)
    
  } else { # LDA or QDA
    # Extract and predict scores
    scoresCal <- scores(object, ncomp = ncomp)
    scoresVal <- predict(object, newdata = newdata, type = "scores", ncomp = ncomp)
    
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


#' @rdname sopls_object
#' @export
R2.sopls <- function(object, estimate, newdata, ncomp = "all", individual = FALSE, ...){
  if (missing(estimate)) {
    ## Select the `best' available estimate
    if (!missing(newdata)) {
      estimate = "test"
    } else {
      if (!is.null(object$validation)) {
        estimate = "CV"
      } else {
        estimate = "train"
      }
    }
  }
  if(!(estimate %in% c("train","CV","test")))
    stop(paste0("'estimate' = ", estimate, " not supported"))
  
  selComp <- pathComp(ncomp, object$decomp$compList)
  
  if(estimate == "CV"){
    if(individual)
      return(object$validation$expl_var_ind[, selComp$hits])
    else
      return(object$validation$expl_var[selComp$hits])
  } else {
    if(estimate == "train"){
      y <- t(object$data$Y)
      yc <- array(object$data$Y, c(dim(object$fitted)))
      sst <- rowSums((y - rowMeans(y))^2)
      ssr <- apply((yc - object$fitted)^2,2:3,sum)
      if(individual){
        r2  <- 1 - ssr/sst
        return(r2[, selComp$hits, drop=FALSE])}
      else{
        r2  <- 1 - colSums(ssr)/sum(sst)
        return(r2[selComp$hits, drop=FALSE])}
    } else {
      # estimate = "test"
      newdata <- model.frame(formula(object), data = newdata)
      preds <- predict(object, newdata, ncomp = object$ncomp, comps = ncomp)
      y <- t(as.matrix(model.response(newdata)))
      yc <- array(t(y), c(dim(preds)))
      sst <- rowSums((y - rowMeans(y))^2)
      ssr <- apply((yc - preds)^2,2:3,sum)
      if(individual){
        r2  <- 1 - ssr/sst
        return(r2)}
      else{
        r2  <- 1 - colSums(ssr)/sum(sst)
        return(r2)}
    }
  }
}


#' @rdname sopls_object
#' @export
RMSEP.sopls <- function(object, estimate, newdata, ncomp = "all", individual = FALSE, ...){
  if (missing(estimate)) {
    ## Select the `best' available estimate
    if (!missing(newdata)) {
      estimate = "test"
    } else {
      if (!is.null(object$validation)) {
        estimate = "CV"
      } else {
        estimate = "train"
      }
    }
  }
  if(!(estimate %in% c("train","CV","test")))
    stop(paste0("'estimate' = ", estimate, " not supported"))
  
  selComp <- pathComp(ncomp, object$decomp$compList)
  
  if(estimate == "CV"){
    if(individual)
      return(object$validation$RMSECV_ind[, selComp$hits])
    else
      return(object$validation$RMSECV[selComp$hits])
  } else {
    if(estimate == "train"){
      y <- t(object$data$Y)
      yc <- array(object$data$Y, c(dim(object$fitted)))
      MSEP <- apply((yc - object$fitted)^2,2:3,mean)
      if(individual){
        rmsep  <- sqrt(MSEP)
        return(rmsep[, selComp$hits, drop=FALSE])}
      else{
        rmsep <- sqrt(colMeans(MSEP))
        return(rmsep[selComp$hits, drop=FALSE])}
    } else {
      # estimate = "test"
      preds <- predict(object, newdata, ncomp = object$ncomp, comps = ncomp)
      newdata <- model.frame(formula(object), data = newdata)
      y <- t(as.matrix(model.response(newdata)))
      yc <- array(t(y), c(dim(preds)))
      MSEP <- apply((yc - preds)^2,2:3,mean)
      if(individual){
        rmsep  <- sqrt(MSEP)
        return(rmsep)}
      else{
        rmsep <- sqrt(colMeans(MSEP))
        return(rmsep)}
    }
  }
}
