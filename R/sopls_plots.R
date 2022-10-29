#' @name sopls_plots
#' @title Scores, loadings and plots for sopls objects
#'
#' @aliases loadings.sopls scores.sopls loadingplot.sopls scoreplot.sopls biplot.sopls
#'
#' @param object \code{sopls} object
#' @param x \code{sopls} object
#' @param comps \code{integer} vector giving components, within block, to plot (see Details regarding combination of blocks).
#' @param ncomp \code{integer} vector giving components from all blocks before \code{block} (see next argument).
#' @param block \code{integer} indicating which block to extract components from.
#' @param y \code{logical} extract Y loadings/scores instead of X loadings/scores (default = FALSE).
#' @param scatter \code{logical} indicating if a scatterplot of loadings should be made (default = TRUE).
#' @param labels \code{character} indicating if "names" or "numbers" should be plot symbols (optional).
#' @param identify \code{logical} for activating \code{identify} to interactively identify points.
#' @param type \code{character} for selecting type of plot to make. Defaults to "p" (points) for scatter plots and "l" (lines) for line plots.
#' @param which \code{character} for selecting type of biplot ("x" = default, "y", "scores", "loadings").
#' @param var.axes \code{logical} indicating if second axes of a biplot should have arrows.
#' @param xlabs \code{character} vector for labelling first set of biplot points (optional).
#' @param ylabs \code{character} vector for labelling second set of biplot points (optional).
#' @param xlab \code{character} text for x labels.
#' @param ylab \code{character} text for y labels.
#' @param lty Vector of line type specifications (see \code{\link{par}} for details).
#' @param lwd \code{numeric} vector of line width specifications.
#' @param pch Vector of point specifications (see \code{\link{points}} for details).
#' @param cex \code{numeric} vector of plot size expansions (see \code{\link{par}} for details).
#' @param col \code{integer} vector of symbol/line colours (see \code{\link{par}} for details).
#' @param legendpos \code{character} indicating legend position (if \code{scatter} is FALSE), e.g. \code{legendpos = "topright"}.
#' @param pretty.xlabels \code{logical} indicating if xlabels should be more nicely plotted (default = TRUE).
#' @param xlim \code{numeric} vector of length two, with the x limits of the plot (optional).
#' @param main \code{character} for setting the main title of a plot.
#' @param plotx \code{locical} or \code{integer}/\code{character}.  Whether to plot the \eqn{X} correlation loadings, optionally which block(s). Defaults to \code{TRUE}.
#' @param ploty \code{logical}.  Whether to plot the \eqn{Y} correlation loadings. Defaults to \code{TRUE}.
#' @param ... further arguments sent to the underlying plot function(s)
#'
#' @description Extraction of \code{scores} and \code{loadings} and adaptation of \code{scoreplot},
#' \code{loadingplot} and \code{biplot} from package \code{pls} for \code{sopls} objects.
#' 
#' @details If \code{comps} is supplied as a \code{list} for \code{scoreplot}, it is assumed that its elements refer to each of the
#' blocks up to block number \code{block}. For instance \code{comps = list(1, 0, 1:2)} will select 1 component from the first
#' block, no components from the second block and the first two components from the last block. This must be 
#' matched by \code{ncomp}, specifying how many components were selected before block number \code{block}.
#'
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results are found in \code{\link{sopls_results}}.
#' 
#' #' @return The score and loading functions return scores and loadings, while plot functions have no return (except use of 'identify').
#' 
#' @examples
#' data(potato)
#' so <- sopls(Sensory ~ Chemical + Compression + NIRraw, data=potato, ncomp=c(5,5,5))
#' 
#' # Loadings
#' loadings(so, ncomp=c(3), block=2)[, 1:3]
#' 
#' # Scores
#' scores(so, block=1)[, 1:4]
#' 
#' # Default plot from first block
#' scoreplot(so)
#' 
#' # Second block with names
#' scoreplot(so, ncomp=c(3), block=2, labels="names")
#' 
#' # Scatterplot matrix
#' scoreplot(so, ncomp=c(3,2), block=3, comps=1:3)
#' 
#' # Combination of blocks (see Details)
#' scoreplot(so, ncomp=c(3,2), block=3, comps=list(1,0,1))
#' 
#' # Default plot from first block
#' loadingplot(so, scatter=TRUE)
#' 
#' # Second block with names
#' loadingplot(so, ncomp=c(3), block=2, labels="names", scatter=TRUE)
#' 
#' # Scatterplot matrix
#' loadingplot(so, ncomp=c(3,2), block=3, comps=1:3, scatter=TRUE)
#' 
#' # Correlation loadings
#' corrplot(so, block=2, ncomp=1)
#' 
#' # Default plot from first block
#' biplot(so)
#' @export
loadings.sopls <- function(object, ncomp = "all", block = 1, y = FALSE, ...){
  if(is.numeric(ncomp) && length(ncomp)!=(block-1))
    stop("Length of 'ncomp' must be one less than 'block'.")
  if(is.character(ncomp))
    selComp <- pathComp(ncomp, object$decomp$compList)
  else
    selComp <- pathComp(c(ncomp,object$max_comps,rep(0,length(object$ncomp)-block)), object$decomp$compList)
  hits <- selComp$hits[!is.na(selComp$hits)]
  blocks  <- object$decomp$changeBlock[hits]
  T <- object$decomp$T[, hits[block==blocks], drop=FALSE]
  X <- scale(object$data$X[[block]], scale=FALSE)
  P <- crossprod(X, T)
  Xcat <- scale(do.call(cbind,object$data$X), scale=FALSE)
  Pcat <- crossprod(Xcat, T)
  xve <- 100 * apply(Pcat^2,2,sum)/sum(Xcat^2)
  if(y)
    P <- object$decomp$Q[, hits[block==blocks],drop=FALSE]
  attr(P, "explvar") <- xve
  xvei <- 100 * apply(P^2,2,sum)/sum(X^2)
  attr(P, "explvar_block") <- xvei
  class(P) <- c("loadings.multiblock","loadings")
  return(P)
}

#' @rdname sopls_plots
#' @export
scores.sopls <- function(object, ncomp = "all", block = 1, y = FALSE, ...){
  if(is.numeric(ncomp) && length(ncomp)!=(block-1))
    stop("Length of 'ncomp' must be one less than 'block'.")
  if(is.character(ncomp))
    selComp <- pathComp(ncomp, object$decomp$compList)
  else
    selComp <- pathComp(c(ncomp,object$max_comps,rep(0,length(object$ncomp)-block)), object$decomp$compList)
  hits <- selComp$hits[!is.na(selComp$hits)]
  blocks  <- object$decomp$changeBlock[hits]
  T <- object$decomp$T[, hits[block==blocks], drop=FALSE]
  X <- scale(object$data$X[[block]], scale=FALSE)
  Xcat <- scale(do.call(cbind,object$data$X), scale=FALSE)
  Pcat <- crossprod(Xcat, T)
  xve <- 100 * apply(Pcat^2,2,sum)/sum(Xcat^2)
  if(y){
    U <- object$data$Y %*% object$decomp$Q[, hits[block==blocks],drop=FALSE]
    if(sum(block==blocks) > 1)
      for(a in 2:sum(block==blocks))
        U[,a] <- U[,a] - T[,1:a] %*% crossprod(T[,1:a], U[,a])
    T <- U
  }
  attr(T, "explvar") <- xve
  class(T) <- c("scores.multiblock","scores")
  T
}

#' @rdname sopls_plots
#' @export
scoreplot.sopls <- function(object, comps = 1:2, ncomp = NULL, block = 1, labels, identify = FALSE,
                            type = "p", xlab, ylab, ...){
  ## Check arguments
  nComps <- length(comps)
  if(nComps == 0) stop("At least one component must be selected.")
  if(length(ncomp) != block-1){
    if(block == 1)
      stop("Wrong length of ncomp. It must specify the number of components\nin all preceding blocks, i.e., ncomp = NULL.")
    else
      if(block == 2)
        stop(paste0("Wrong length of ncomp. It must specify the number of components\nin all preceding blocks, e.g., ncomp = ", object$ncomp[1], ")."))
      else
        stop(paste0("Wrong length of ncomp. It must specify the number of components\nin all preceding blocks, e.g., ncomp = c(",paste(object$ncomp[1:(block-1)], collapse=","),")."))
  }
  ## Get the scores
  if (is.matrix(object)) {
    ## Assume this is already a score matrix
    S <- object
  } else {
    ## Try to get the scores
    if(is.list(comps)){
      if(length(comps) != block)
        stop("When supplying 'comps' as a list to plot a combination of components\nfrom different blocks, 'block' must equal length of 'comps'.")
      S <- matrix(0,nrow=nrow(object$data$Y),0)
      evar <- numeric()
      for(i in 1:length(comps)){
        if(!all(comps[[i]]==0))
          if(i>1)
            s <- scores(object, ncomp = ncomp[1:(i-1)], block = i)
          else
            s <- scores(object, ncomp = NULL, block = i)
          evar <- cbind(evar, attr(s,'explvar')[comps[[i]]])
          if(dim(s)[2]==1 && length(comps[[i]])>1){
            comps[[i]] <- comps[[i]][1]
            warning(paste0("Only one component in block ",i," but multiple components selected"))
          }
          nComps <- sum(unlist(comps)!=0)
          S <- cbind(S,s[,comps[[i]], drop = FALSE])
      }
    } else
      S <- scores(object, ncomp = ncomp, block = block)
    if (is.null(S))
      stop("`", deparse(substitute(object)), "' has no scores.")
  }
  if(dim(S)[2]==1 && length(comps)>1){
    comps <- comps[1]
    nComps <- 1
    warning(paste0("Only one component in block but multiple components selected: comps = c(", paste(comps,collapse=","), ")"))
  }
  if(!is.list(comps)){
    evar <- attr(S,'explvar')[comps]
    S <- S[,comps, drop = FALSE]
  }
  varlab <- paste(colnames(S), " (", format(evar, digits = 2, trim = TRUE),
                  " %)", sep = "")
  if (!missing(labels)) {
    ## Set up point labels
    if (length(labels) == 1) {
      labels <- switch(match.arg(labels, c("names", "numbers")),
                       names = rownames(S),
                       numbers = 1:nrow(S)
      )
    }
    labels <- as.character(labels)
    type <- "n"
  }
  if (nComps <= 2) {
    if (nComps == 1) {
      ## One component versus index
      if (missing(xlab)) xlab <- "observation"
      if (missing(ylab)) ylab <- varlab
    } else {
      ## Second component versus first
      if (missing(xlab)) xlab <- varlab[1]
      if (missing(ylab)) ylab <- varlab[2]
    }
    plot(S, xlab = xlab, ylab = ylab, type = type, ...)
    if (!missing(labels)) text(S, labels, ...)
    if (isTRUE(identify)) {
      if (!is.null(rownames(S))) {
        identify(S, labels = rownames(S))
      } else {
        identify(S)
      }
    }
  } else {
    ## Pairwise scatterplots of several components
    panel <- if (missing(labels))
      function(x, y, ...) points(x, y, type = type, ...) else
        function(x, y, ...) text(x, y, labels = labels, ...)
    pairs(S, labels = varlab, panel = panel, ...)
  }
}

#' @rdname sopls_plots
#' @export
loadingplot.sopls <- function(object, comps = 1:2, ncomp = NULL, block = 1, scatter = TRUE, labels,
                              identify = FALSE, type, lty, lwd = NULL, pch,
                              cex = NULL, col, legendpos, xlab, ylab,
                              pretty.xlabels = TRUE, xlim, ...)
{
  ## Check arguments
  nComps <- length(comps)
  if (nComps == 0) stop("At least one component must be selected.")
  if(length(ncomp) != block-1){
    if(block == 1)
      stop("Wrong length of ncomp. It must specify the number of components\nin all preceding blocks, i.e., ncomp = NULL.")
    else
      if(block == 2)
        stop(paste0("Wrong length of ncomp. It must specify the number of components\nin all preceding blocks, e.g., ncomp = ", object$ncomp[1], ")."))
    else
      stop(paste0("Wrong length of ncomp. It must specify the number of components\nin all preceding blocks, e.g., ncomp = c(",paste(object$ncomp[1:(block-1)], collapse=","),")."))
  }
  if (!missing(type) &&
      (length(type) != 1 || is.na(nchar(type, "c")) || nchar(type, "c") != 1))
    stop("Invalid plot type.")
  ## Get the loadings
  if (is.matrix(object)) {
    ## Assume this is already a loading matrix
    L <- object
  } else {
    ## Try to get the loadings:
    L <- loadings(object, ncomp = ncomp, block = block)
    if (is.null(L))
      stop("`", deparse(substitute(object)), "' has no loadings.")
  }
  if(dim(L)[2]==1 && length(comps)>1){
    comps <- comps[1]
    nComps <- 1
    warning(paste0("Only one component in block but multiple components selected: comps = c(", paste(comps,collapse=","), ")"))
  }
  evar <- attr(L,'explvar')[comps]
  L <- L[,comps, drop = FALSE]
  varlab <- paste(colnames(L), " (", format(evar, digits = 2, trim = TRUE),
                  " %)", sep = "")
  if (isTRUE(scatter)) {
    ## Scatter plots
    if (missing(type)) type <- "p"
    if (!missing(labels)) {
      ## Set up point/tick mark labels
      if (length(labels) == 1) {
        labels <- switch(match.arg(labels, c("names", "numbers")),
                         names = {
                           if (is.null(rnames <- rownames(L))) {
                             stop("The loadings have no row names.")
                           } else {
                             rnames
                           }},
                         numbers = 1:nrow(L)
        )
      }
      labels <- as.character(labels)
      type <- "n"
    }
    if (missing(lty)) lty <- NULL
    if (missing(pch)) pch <- NULL
    if (missing(col)) col <- par("col") # `NULL' means `no colour'
    if (nComps <= 2) {
      if (nComps == 1) {
        ## One component versus index
        if (missing(xlab)) xlab <- "variable"
        if (missing(ylab)) ylab <- varlab
      } else {
        ## Second component versus first
        if (missing(xlab)) xlab <- varlab[1]
        if (missing(ylab)) ylab <- varlab[2]
      }
      plot(L, xlab = xlab, ylab = ylab, type = type, lty = lty,
           lwd = lwd, pch = pch, cex = cex, col = col, ...)
      if (!missing(labels)) text(L, labels, cex = cex, col = col, ...)
      if (isTRUE(identify))
        identify(L, labels = paste(1:nrow(L), rownames(L), sep = ": "))
    } else {
      ## Pairwise scatterplots of several components
      panel <- if (missing(labels)) {
        function(x, y, ...)
          points(x, y, type = type, lty = lty, lwd = lwd,
                 pch = pch, col = col, ...)
      } else {
        function(x, y, ...)
          text(x, y, labels = labels, col = col, ...)
      }
      pairs(L, labels = varlab, panel = panel, cex = cex, ...)
    }
  } else {                            # if (isTRUE(scatter))
    ## Line plots
    if (missing(type)) type <- "l"
    if (missing(lty))  lty  <- 1:nComps
    if (missing(pch))  pch  <- 1:nComps
    if (missing(col))  col  <- 1:nComps
    if (missing(xlab)) xlab <- "variable"
    if (missing(ylab)) ylab <- "loading value"
    xnum <- 1:nrow(L)
    if (missing(labels)) {
      xaxt <- par("xaxt")
    } else {
      xaxt <- "n"
      if (length(labels) == 1) {
        xnam <- rownames(L)
        switch(match.arg(labels, c("names", "numbers")),
               names = {        # Simply use the names as is
                 labels <- xnam
               },
               numbers = {      # Try to use them as numbers
                 if (length(grep("^[-0-9.]+[^0-9]*$", xnam)) ==
                     length(xnam)) {
                   ## Labels are on "num+text" format
                   labels <- sub("[^0-9]*$", "", xnam)
                   if (isTRUE(pretty.xlabels)) {
                     xnum <- as.numeric(labels)
                     xaxt <- par("xaxt")
                   }
                 } else {
                   stop("Could not convert variable names to numbers.")
                 }
               }
        )
      } else {
        labels <- as.character(labels)
      }
    }
    if (missing(xlim)) xlim <- xnum[c(1, length(xnum))] # Needed for reverted scales
    matplot(xnum, L, xlab = xlab, ylab = ylab, type = type,
            lty = lty, lwd = lwd, pch = pch, cex = cex, col = col,
            xaxt = xaxt, xlim = xlim, ...)
    if (!missing(labels) && xaxt == "n") {
      if (isTRUE(pretty.xlabels)) {
        ticks <- axTicks(1)
        ticks <- ticks[ticks >= 1 & ticks <= length(labels)]
      } else {
        ticks <- 1:length(labels)
      }
      axis(1, ticks, labels[ticks], ...)
    }
    if (!missing(legendpos)) {
      ## Are we plotting lines?
      dolines <- type %in% c("l", "b", "c", "o", "s", "S", "h")
      ## Are we plotting points?
      dopoints <- type %in% c("p", "b", "o")
      if (length(lty) > nComps) lty <- lty[1:nComps]
      do.call("legend", c(list(legendpos, varlab, col = col),
                          if (dolines) list(lty = lty, lwd = lwd),
                          if (dopoints) list(pch = pch, pt.cex = cex,
                                             pt.lwd = lwd)))
    }
    if (isTRUE(identify))
      identify(c(row(L)), c(L),
               labels = paste(c(col(L)), rownames(L), sep = ": "))
  }                                   # if (isTRUE(scatter))
}

#' @rdname sopls_plots
#' @export
#' 
corrplot.sopls <- function(object, comps=1:2, ncomp = NULL, block = 1, labels=TRUE, col=1:5, 
                           plotx=TRUE, ploty=TRUE, ...){
  ## Check arguments
  nComps <- length(comps)
  if(nComps == 0) stop("At least one component must be selected.")
  if(length(ncomp) != block-1){
    if(block == 1)
      stop("Wrong length of ncomp. It must specify the number of components\nin all preceding blocks, i.e., ncomp = NULL.")
    else
      if(block == 2)
        stop(paste0("Wrong length of ncomp. It must specify the number of components\nin all preceding blocks, e.g., ncomp = ", object$ncomp[1], ")."))
    else
      stop(paste0("Wrong length of ncomp. It must specify the number of components\nin all preceding blocks, e.g., ncomp = c(",paste(object$ncomp[1:(block-1)], collapse=","),")."))
  }
  ## Get the scores
  if (is.matrix(object)) {
    ## Assume this is already a score matrix
    S <- object
  } else {
    ## Try to get the scores
    S <- scores(object, ncomp = ncomp, block = block)
    if (is.null(S))
      stop("`", deparse(substitute(object)), "' has no scores.")
  }
  if(dim(S)[2]==1 && length(comps)>1){
    comps <- comps[1]
    nComps <- 1
    warning(paste0("Only one component in block but multiple components selected: comps = c(", paste(comps,collapse=","), ")"))
  }
  pls::corrplot(S[0,], plotx=FALSE, ploty=FALSE, comps=comps, ...)
  if(is.logical(plotx) && plotx){
    plotx <- 1:length(object$data$X)
  }
  if(!(is.logical(plotx) && !plotx)){
    for(i in plotx){
      if(labels){
        text(cor(object$data$X[[i]], S[,comps]), labels=colnames(object$data$X[[i]]), col=col[i])
      } else {
        points(cor(object$data$X[[i]], S[,comps]), col=col[i])
      }
    }
  }
  if(ploty && !is.null(object$data$Y)){
    if(labels)
      text(cor(object$data$Y, S[,comps]), labels=colnames(object$data$Y), col=col[length(object$data$X)+1])
    else
      points(cor(object$data$Y, S[,comps]), col=col[length(object$data$X)+1])
  }
}


#' @export
#' @rdname sopls_plots
biplot.sopls <- function(x, comps = 1:2, ncomp = "all", block = 1, which = c("x", "y", "scores", "loadings"),
                              var.axes = FALSE, xlabs, ylabs, main, ...)
{
  if (length(comps) != 2) stop("Exactly 2 components must be selected.")
  which <- match.arg(which)
  switch(which,
         x = {
           objects <- scores(x, ncomp = ncomp, block = block)
           vars <- loadings(x, ncomp = ncomp, block = block)
           title <- "X scores and X loadings"
         },
         y = {
           objects <- scores(x, ncomp = ncomp, block = block, y = TRUE)
           vars <- loadings(x, ncomp = ncomp, block = block, y = TRUE)
           title <- "Y scores and Y loadings"
         },
         scores = {
           objects <- scores(x, ncomp = ncomp, block = block)
           vars <- scores(x, ncomp = ncomp, block = block, y = TRUE)
           title <- "X scores and Y scores"
         },
         loadings = {
           objects <- loadings(x, ncomp = ncomp, block = block)
           vars <- loadings(x, ncomp = ncomp, block = block, y = TRUE)
           title <- "X loadings and Y loadings"
         }
  )
  if (is.null(objects) || is.null(vars))
    stop("'x' lacks the required scores/loadings.")
  ## Build a call to `biplot'
  mc <- match.call()
  mc$comps <- mc$which <- NULL
  mc$x <- objects[,comps, drop = FALSE]
  mc$y <- vars[,comps, drop = FALSE]
  mc$block <- mc$ncomp <- NULL
  if (missing(main)) mc$main <- title
  if (missing(var.axes)) mc$var.axes = FALSE
  if (!missing(xlabs) && isFALSE(xlabs))
    mc$xlabs <- rep("o", nrow(objects))
  if (!missing(ylabs) && isFALSE(ylabs))
    mc$ylabs <- rep("o", nrow(vars))
  mc[[1]] <- as.name("biplot")
  ## Evaluate the call:
  eval(mc, parent.frame())
}