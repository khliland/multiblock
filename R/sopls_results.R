#' @name sopls_results
#' @title Result functions for SO-PLS models
#'
#' @aliases predict.sopls coef.sopls print.sopls summary.sopls pcp.sopls R2.sopls classify.sopls RMSEP.sopls cvanova.sopls cvanova
#' @param object A \code{sopls} object.
#' @param x A \code{sopls} object.
#' @param newdata Optional new data with the same types of predictor blocks as the ones used for fitting the object.
#' @param ncomp An \code{integer} vector giving the exact components to apply.
#' @param type A \code{character} for \code{predict} indicating if responses or scores should be predicted (default = "response", or "scores"), for \code{summary} indicating which type of explained variance to compute (default = "train", alternative = "CV").
#' @param na.action Function determining what to do with missing values in \code{newdata}.
#' @param intercept A \code{logical} indicating if coefficients for the intercept should be included (default = FALSE).
#' @param what A \code{character} indicating if summary should include all, validation or training.
#' @param digits The number of digits used for printing.
#' @param print.gap Gap between columns when printing.
#' @param classes A \code{character} vector of class labels.
#' @param LQ A \code{character} indicating if 'max' (maximum score value), 'lda' or 'qda' should be used when classifying.
#' @param estimate A \code{character} indicating if 'train', 'CV' or 'test' results should be displayed.
#' @param individual A \code{logical} indicating if results for individual responses should be displayed.
#' @param X A \code{list} of data blocks.
#' @param pred An object holding the CV-predicted values (\code{sopls}, \code{matrix} or \code{list} of vectors)
#' @param true A \code{numeric} of true response values for CVANOVA.
#' @param absRes A \code{logical} indicating if absolute (TRUE) or squared (FALSE) residuals should be computed.
#' @param comps An \code{integer} vector giving the exact components to apply.
#' @param ... Additional arguments. Currently not implemented.
#'
#' @return Returns depend on method used, e.g. \code{predict.sopls} returns predicted responses 
#' or scores depending on inputs, \code{coef.sopls} return regression coefficients, while print and summary methods return the object invisibly.
#' 
#' @description Standard result functions for SO-PLS (\code{\link{sopls}}).
#' 
#' @details The parameter \code{ncomp} controls
#' which components to apply/extract, resulting in the sequence of components leading up to the specific choice, i.e.
#' \code{ncomp = c(2,2,1)} results in the sequence 1,0,0; 2,0,0; 2,1,0; 2,2,0; 2,2,1.
#' Usage of the functions are shown using generics in the examples below. 
#' Prediction, regression coefficients, object printing and summary are available through: 
#' \code{predict.sopls}, \code{coef.sopls}, \code{print.sopls} and \code{summary.sopls}.
#' Explained variances and RMSEP are available through \code{R2.sopls} and \code{RMSEP.sopls}.
#' Principal components of predictions are available through \code{pcp.sopls}. Finally, there is work in progress on classifcation
#' support through \code{classify.sopls}.
#' 
#' @references Jørgensen K, Mevik BH, Næs T. Combining designed experiments with several blocks of spectroscopic data. Chemometr Intell Lab Syst. 2007;88(2): 154–166.
#'
#' @examples
#' data(potato)
#' mod <- sopls(Sensory[,1] ~ ., data = potato[c(1:3,9)], ncomp = 5, subset = 1:20)
#' testset <- potato[-(1:20),]; testset$Sensory <- testset$Sensory[,1,drop=FALSE]
#' predict(mod, testset, ncomp=c(2,1,2))
#' dim(coef(mod, ncomp=c(3,0,1))) # <variables x responses x components>
#' R2(mod, ncomp = c(4,1,2))
#' print(mod)
#' summary(mod)
#' 
#' # PCP from sopls object
#' modMulti <- sopls(Sensory ~ ., data = potato[c(1:3,9)], ncomp = 5, validation = "CV", segment = 5)
#' (PCP <- pcp(modMulti, c(2,1,2)))
#' scoreplot(PCP)
#' 
#' # PCP from matrices
#' preds <- modMulti$validation$Ypred[,,"2,1,2"]
#' PCP_default <- pcp(preds, potato[1:3])
#' 
#' # CVANOVA
#' modCV <- sopls(Sensory[,1] ~ ., data = potato[c(1:3,9)], ncomp = 5, validation = "CV", segment = 5)
#' summary(cva <- cvanova(modCV, "2,1,2"))
#' plot(cva)
#' 
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for plotting are found in \code{\link{sopls_plots}}.
#' @importFrom mixlm simple.glht cld
#' @importFrom stats anova lm
#' @export
predict.sopls <- function(object, newdata, ncomp = object$ncomp,
                          type = c("response", "scores"), na.action = na.pass, ...){
  if (missing(newdata) || is.null(newdata)){
    newdata <- object$data$X
  } else {
    newdata <- model.frame(formula(object), data = newdata)
    newdata <- newdata[-1]
  }
  
  type <- match.arg(type)
  
  if(!is.character(ncomp) && sum(ncomp) > object$max_comp)
    stop(paste0("Selected components (",paste(ncomp,collapse=","),") is outside range of fitted model."))
  if (type == "response") {
    return(sopls_prediction(object, newdata, ncomp, FALSE))
  } else {
    return(sopls_prediction(object, newdata, ncomp, TRUE))
  }
}

#' @rdname sopls_results
#' @export
coef.sopls <- function(object, ncomp = object$ncomp, intercept = FALSE,
                       ...)
{
  if(sum(ncomp) > object$max_comp)
    stop(paste0("Selected components (",paste(ncomp,collapse=","),") is outside range of fitted model."))
  
  X <- object$data$X
  Y <- object$data$Y
  n    <- dim(X[[1]])[1]
  
  nblock  <- length(X)
  nresp   <- dim(Y)[2]
  selComp <- pathComp(ncomp, object$decomp$compList)
  tot_comp <- nrow(selComp$path)
  
  Cr <- Crval <- 0
  for(i in 1:nblock){
    if(ncomp[i]>0){
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

#' @rdname sopls_results
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

#' @rdname sopls_results
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


#' @rdname sopls_results
#' @export
classify <- function(object, ...) UseMethod("classify")

#' @rdname sopls_results
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


#' @rdname sopls_results
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
      preds <- predict(object, newdata, ncomp = ncomp)
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


#' @rdname sopls_results
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
      preds <- predict(object, newdata, ncomp = ncomp)
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

#' @rdname sopls_results
#' @export
pcp <- function (object, ...) {
  UseMethod("pcp", object)
}

#' @rdname sopls_results
#' @export
pcp.sopls <- function(object, ncomp, ...){
  if(is.null(object$validation))
    stop("The 'sopls' object has been created without validation")
  if(missing(ncomp))
    stop("The argument 'ncomp' must be supplied.")

  compName <- paste0(ncomp, collapse = ",")
  preds <- object$validation$Ypred[,,compName,drop=FALSE]
  dim(preds) <- dim(object$validation$Ypred)[1:2]
  dimnames(preds) <- dimnames(object$validation$Ypred)[1:2]
  PCP <- pca(preds, ncomp = min(ncol(preds),nrow(preds)-1))
  PCP$data = object$data
  PCP$loadings <- PCP$loadings * rep(sqrt(colSums(PCP$scores^2)),each=nrow(PCP$loadings))
  PCP$scores <- PCP$scores / rep(sqrt(colSums(PCP$scores^2)),each=nrow(PCP$scores))
  PCP$blockLoadings <- lapply(lapply(object$data$X, function(x)x-rep(colMeans(x),each=nrow(x))),function(x)crossprod(x,PCP$scores))
  PCP$info <- list(method = "Principal Components of Predictions", 
                   scores = "Scores", loadings = "Loadings",
                   blockScores = "not used", blockLoadings = "Block loadings")
  PCP$call <- match.call()
  class(PCP) <- c('multiblock','list')
  return(PCP)
  
}

#' @rdname sopls_results
#' @export
pcp.default <- function(object, X, ...){
  # Assumes object to be a matrix/data.frame of predictions (possibly CV)
  PCP <- pca(object, ncomp = min(ncol(object),nrow(object)-1))
  PCP$loadings <- PCP$loadings * rep(sqrt(colSums(PCP$scores^2)),each=nrow(PCP$loadings))
  PCP$scores <- PCP$scores / rep(sqrt(colSums(PCP$scores^2)),each=nrow(PCP$scores))
  PCP$info <- list(method = "Principal Components of Predictions", 
                   scores = "Scores", loadings = "Loadings",
                   blockScores = "not used", blockLoadings = "not used")
  if(!missing(X)){
    PCP$blockLoadings <- lapply(lapply(X, function(x)x-rep(colMeans(x),each=nrow(x))),function(x)crossprod(x,PCP$scores))
    PCP$info$blockLoadings <- "Block loadings"
    PCP$data <- X
  }
  PCP$call <- match.call()
  class(PCP) <- c('multiblock','list')
  return(PCP)
}

#' @rdname sopls_results
#' @export
cvanova <- function (pred, ...) {
  UseMethod("cvanova", pred)
}

#' @rdname sopls_results
#' @export
cvanova.default <- function(pred, true, absRes = TRUE, ...){
  if(!is.null(dim(true)) && dim(true)[2]>1)
    stop("'cvanova' only supports single response models")
  if(inherits(pred, 'list'))
    pred <- do.call(cbind, pred)
  modNames <- 1:ncol(pred)
  if(!is.null(colnames(pred)))
    modNames <- colnames(pred)
  DF <- data.frame('Model' = factor(rep(modNames, each=length(true))), 'Object' = factor(1:length(true)), 'Residual' = c(abs(true-pred)))
  if(!absRes)
    DF[['Residual']] <- DF[['Residual']]^2
  mod <- list()
  mod$lm <- lm(Residual ~ Model + Object, data = DF)
  mod$anova <- anova(mod$lm)
  mod$tukey <- simple.glht(mod$lm,'Model')
  class(mod) <- c('cvanova','list')
  return(mod)
}

#' @rdname sopls_results
#' @export
cvanova.sopls <- function(pred, comps, absRes = TRUE, ...){
  if(is.null(pred$validation))
    stop("The 'sopls' object has been created without validation")
  if(!is.null(dim(pred$data$Y)) && dim(pred$data$Y)[2]>1)
    stop("'cvanova' only supports single response models")
  compSplit <- strsplit(comps, ",")[[1]]
  nBlock <- length(compSplit)
  compVec <- character(nBlock)
  for(i in 1:nBlock){
    zeros <- rep("0", nBlock-i)
    compVec[i] <- paste0(c(compSplit[1:i],zeros), collapse = ",")
  }
  preds <- pred$validation$Ypred[,1,compVec]
  trues <- c(pred$data$Y)
  return(cvanova(preds, trues, absRes))
}

#' @rdname sopls_results
#' @export
print.cvanova <- function(x, ...){
  print(x$lm)
}

#' @rdname sopls_results
#' @export
summary.cvanova <- function(object, ...){
  print(object$anova)
  print(cld(object$tukey))
}

#' @rdname sopls_results
#' @export
plot.cvanova <- function(x, ...){
  plot(x$tukey)
}

#' @rdname sopls_results
#' @export
residuals.sopls <- function(object, ...){
  res <- object$fitted
  for(j in 1:dim(res)[3])
    res[,,j] <- object$data$Y - object$fitted[,,j]
  return(res)
}
