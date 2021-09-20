#' @name maage
#' @title Måge plot
#'
#' @param object An SO-PLS model (\code{sopls} object)
#' @param compSeq Integer vector giving the sequence of previous components chosen for \code{maageSeq} (see example).
#' @param expl_var Logical indicating if explained variance (default) or RMSECV should be displayed.
#' @param pure.trace Logical indicating if single block solutions should be traced in the plot.
#' @param pch Scalar or symbol giving plot symbol.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param xlim Plot limits for x-axis (numeric vector).
#' @param ylim Plot limits for y-axis (numeric vector).
#' @param cex.text Text scaling (scalar) for better readability of plots.
#' @param col Line colour in plot.
#' @param col.block Line colours for blocks (default = c('red','blue','darkgreen','purple','black'))
#' @param ... Additional arguments to \code{plot}.
#' 
#' @description Måge plot for SO-PLS (\code{\link{sopls}}) cross-validation visualisation. 
#' 
#' @details This function can either be used 
#' for global optimisation across blocks or sequential optimisation, using \code{maageSeq}.
#' The examples below show typical usage.
#' 
#' @return The \code{maage} plot has no return.
#'
#' @examples
#' data(wine)
#' ncomp <- unlist(lapply(wine, ncol))[-5]
#' so.wine <- sopls(`Global quality` ~ ., data=wine, ncomp=ncomp, 
#'             max_comps=10, validation="CV", segments=10)
#' maage(so.wine)
#' 
#' # Sequential search for optimal number of components per block
#' old.par <- par(mfrow=c(2,2), mar=c(3,3,0.5,1), mgp=c(2,0.7,0))
#' maageSeq(so.wine)
#' maageSeq(so.wine, 2)
#' maageSeq(so.wine, c(2,1))
#' maageSeq(so.wine, c(2,1,1))
#' par(old.par)
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' @export
maage <- function(object, expl_var=TRUE, pure.trace=FALSE, pch=20, xlab='# components', 
                  ylab=ifelse(expl_var,'Explained variance (%)','RMSECV'), 
                  xlim=NULL, ylim=NULL, cex.text=0.8, ...){
  x <- rowSums(object$decomp$compList)
  if(is.null(xlim)){
    xlim <- c(0,max(x)+0.7)
  }
  nblock <- dim(object$decomp$compList)[2]
  if(expl_var){ # Explained variance
    if(is.null(ylim)){
      ylim <- c(min(object$validation$expl_var*100),100)
    }
    plot(x, object$validation$expl_var*100, pch=pch, 
         xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
    text(rowSums(object$decomp$compList), object$validation$expl_var*100, dimnames(object$validation$Ypred)[[3]], pos=4, cex=cex.text, ...)
    if(pure.trace){ # Trace of single block
      for(i in 1:nblock){
        inds <- which(rowSums(object$decomp$compList[,-i,drop=FALSE])==0)
        lines(object$decomp$compList[inds,i], object$validation$expl_var[inds]*100, lty=i)
      }
    }
  } else { # RMSECV
    if(is.null(ylim)){
      ylim <- range(object$validation$RMSECV)
    }
    plot(x, object$validation$RMSECV, pch=pch, 
         xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
    text(rowSums(object$decomp$compList), object$validation$RMSECV, dimnames(object$validation$Ypred)[[3]], pos=4, cex=cex.text, ...)
    if(pure.trace){ # Trace of single block
      for(i in 1:nblock){
        inds <- which(rowSums(object$decomp$compList[,-i,drop=FALSE])==0)
        lines(object$decomp$compList[inds,i], object$validation$RMSECV[inds], lty=i)
      }
    }
  }
}

#' @rdname maage
#' @export
maageSeq <- function(object, compSeq=TRUE, expl_var=TRUE, pch=20, xlab='# components', 
                     ylab=ifelse(expl_var,'Explained variance (%)','RMSECV'), 
                     xlim=NULL, ylim=NULL, cex.text=0.8, col='gray', col.block=c('red','blue','darkgreen','purple','black','red','blue','darkgreen'), ...){
  x <- rowSums(object$decomp$compList)
  if(is.null(xlim)){
    xlim <- c(0,max(x)+0.7)
  }
  if(!is.logical(compSeq) && sum(compSeq)>object$max_comps)
    stop('Selected components outside of component range')
  if(length(compSeq)>length(object$ncomp))
    stop('Too many blocks specified')
  nblock <- dim(object$decomp$compList)[2]
  if(expl_var){ # Explained variance
    if(is.null(ylim)){
      ylim <- c(min(object$validation$expl_var*100),100)
    }
    plot(x, object$validation$expl_var*100, pch=pch, 
         xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=col, ...)
    text(rowSums(object$decomp$compList), object$validation$expl_var*100, dimnames(object$validation$Ypred)[[3]], pos=4, cex=cex.text, col=col, ...)
    
    # Trace of Block 1
    inds <- which(rowSums(object$decomp$compList[,-1,drop=FALSE])==0)
    lines(object$decomp$compList[inds,1], object$validation$expl_var[inds]*100, lty=1, col=col.block[1])
    points(object$decomp$compList[inds,1], object$validation$expl_var[inds]*100, lty=1, col=col.block[1],pch=pch)
    text(object$decomp$compList[inds,1], object$validation$expl_var[inds]*100, dimnames(object$validation$Ypred)[[3]][inds], pos=4, cex=cex.text, col=col.block[1], ...)
    if(is.numeric(compSeq)){ 
      traceLimit <- !logical(dim(object$decomp$compList)[1])
      # Trace of next blocks
      for(i in 2:(length(compSeq)+1)){
        traceLimit <- traceLimit & object$decomp$compList[,i-1] == compSeq[i-1]
        inds <- which(rowSums(object$decomp$compList[,-(1:i),drop=FALSE])==0 & traceLimit)
        if(i<=dim(object$decomp$compList)[2]){
          lines(object$decomp$compList[inds,i]+sum(compSeq[1:(i-1)]), object$validation$expl_var[inds]*100, lty=1, col=col.block[i])
          points(object$decomp$compList[inds,i]+sum(compSeq[1:(i-1)]), object$validation$expl_var[inds]*100, lty=1, col=col.block[i], pch=pch)
          text(object$decomp$compList[inds,i]+sum(compSeq[1:(i-1)]), object$validation$expl_var[inds]*100, dimnames(object$validation$Ypred)[[3]][inds], pos=4, cex=cex.text, col=col.block[i], ...)
        } else {
          points(sum(compSeq[1:(i-1)]), object$validation$expl_var[inds]*100, lty=1, col=col.block[i], pch=pch)
          text(sum(compSeq[1:(i-1)]), object$validation$expl_var[inds]*100, dimnames(object$validation$Ypred)[[3]][inds], pos=4, cex=cex.text, col=col.block[i], ...)
        }
      }
    }
  } else { # RMSECV
    if(is.null(ylim)){
      ylim <- range(object$validation$RMSECV)
    }
    plot(x, object$validation$RMSECV, pch=pch, 
         xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=col, ...)
    text(rowSums(object$decomp$compList), object$validation$RMSECV, dimnames(object$validation$Ypred)[[3]], pos=4, cex=cex.text, col=col, ...)
    # Trace of Block 1
    inds <- which(rowSums(object$decomp$compList[,-1,drop=FALSE])==0)
    lines(object$decomp$compList[inds,1], object$validation$RMSECV[inds], lty=1, col=col.block[1])
    points(object$decomp$compList[inds,1], object$validation$RMSECV[inds], lty=1, col=col.block[1], pch=pch)
    text(object$decomp$compList[inds,1], object$validation$RMSECV[inds], dimnames(object$validation$Ypred)[[3]][inds], pos=4, cex=cex.text, col=col.block[1], ...)
    if(is.numeric(compSeq)){ 
      traceLimit <- !logical(dim(object$decomp$compList)[1])
      # Trace of next blocks
      for(i in 2:(length(compSeq)+1)){
        traceLimit <- traceLimit & object$decomp$compList[,i-1] == compSeq[i-1]
        inds <- which(rowSums(object$decomp$compList[,-(1:i),drop=FALSE])==0 & traceLimit)
        if(i<=dim(object$decomp$compList)[2]){
          lines(object$decomp$compList[inds,i]+sum(compSeq[1:(i-1)]), object$validation$RMSECV[inds], lty=1, col=col.block[i])
          points(object$decomp$compList[inds,i]+sum(compSeq[1:(i-1)]), object$validation$RMSECV[inds], lty=1, col=col.block[i],pch=pch)
          text(object$decomp$compList[inds,i]+sum(compSeq[1:(i-1)]), object$validation$RMSECV[inds], dimnames(object$validation$Ypred)[[3]][inds], pos=4, cex=cex.text, col=col.block[i], ...)
        } else {
          points(sum(compSeq[1:(i-1)]), object$validation$RMSECV[inds], lty=1, col=col.block[i],pch=pch)
          text(sum(compSeq[1:(i-1)]), object$validation$RMSECV[inds], dimnames(object$validation$Ypred)[[3]][inds], pos=4, cex=cex.text, col=col.block[i], ...)
        }
      }
    }
  }
}
