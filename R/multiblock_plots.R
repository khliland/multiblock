#' @name multiblock_plots
#' @title Plot Functions for Multiblock Objects
#' @aliases scoreplot.multiblock loadingplot.multiblock loadingweightplot
#' 
#' @description Plotting procedures for \code{multiblock} objects.
#' 
#' @details Plot functions for \code{scores}, \code{loadings} and \code{loading.weights} based
#' on the functions found in the \code{pls} package.
#' 
#' @param object \code{multiblock} object.
#' @param x \code{multiblock} object.
#' @param comps \code{integer} vector giving components, within block, to plot.
#' @param block \code{integer/character} for block selection.
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
#' @param main \code{character} text for main header.
#' @param lty Vector of line type specifications (see \code{\link{par}} for details).
#' @param lwd \code{numeric} vector of line width specifications.
#' @param pch Vector of point specifications (see \code{\link{points}} for details).
#' @param cex \code{numeric} vector of plot size expansions (see \code{\link{par}} for details).
#' @param col \code{integer} vector of symbol/line colours (see \code{\link{par}} for details).
#' @param legendpos \code{character} indicating legend position (if \code{scatter} is FALSE), e.g. \code{legendpos = "topright"}.
#' @param pretty.xlabels \code{logical} indicating if xlabels should be more nicely plotted (default = TRUE).
#' @param xlim \code{numeric} vector of length two, with the x limits of the plot (optional).
#' @param ... Not implemented.
#' 
#' @return These plotting routines only generate plots and return no values.
#' 
#' @examples 
#' data(wine)
#' sc <- sca(wine[c('Smell at rest', 'View', 'Smell after shaking')], ncomp = 4)
#' plot(loadings(sc, block = 1), labels = "names", scatter = TRUE)
#' 
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results are found in \code{\link{multiblock_results}}.
#' @importFrom graphics axTicks matplot pairs
#' @export
scoreplot.multiblock <- function(object, comps = 1:2, block = 0, labels, identify = FALSE,
                                 type = "p", xlab, ylab, main, ...){
  ## Check arguments
  nComps <- length(comps)
  if (nComps == 0) stop("At least one component must be selected.")
  ## Get the scores
  if (is.matrix(object)) {
    ## Assume this is already a score matrix
    S <- object
  } else {
    ## Try to get the scores
    S <- scores(object, block = block)
    if(block != 0){
      scoreType <- object$info$blockScores
      if(!is.null(bn <- names(object$blockScores)))
        if(is.numeric(block))
          scoreType <- paste0(scoreType, ", ", bn[block])
      else
        scoreType <- paste0(scoreType, ", ", block)
    }
    else
      scoreType <- object$info$scores
    if (is.null(S))
      stop("`", deparse(substitute(object)), "' has no scores.")
  }
  evar <- attr(S,'explvar')[comps]
  S <- S[,comps, drop = FALSE]
  if(is.null(evar))
    varlab <- colnames(S)
  else
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
  if(missing(main))
    main <- ifelse(scoreType != "Not used", scoreType, "")
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
    plot(S, xlab = xlab, ylab = ylab, type = type, main = main, ...)
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
    pairs(S, labels = varlab, panel = panel, main = main, ...)
  }
}

#' @rdname multiblock_plots 
#' @export
loadingplot.multiblock <- function(object, comps = 1:2, block = 0, scatter = TRUE, labels,
                                   identify = FALSE, type, lty, lwd = NULL, pch,
                                   cex = NULL, col, legendpos, xlab, ylab, main,
                                   pretty.xlabels = TRUE, xlim, ...)
{
  ## Check arguments
  nComps <- length(comps)
  if (nComps == 0) stop("At least one component must be selected.")
  if (!missing(type) &&
      (length(type) != 1 || is.na(nchar(type, "c")) || nchar(type, "c") != 1))
    stop("Invalid plot type.")
  ## Get the loadings
  if (is.matrix(object)) {
    ## Assume this is already a loading matrix
    L <- object
  } else {
    ## Try to get the loadings:
    L <- loadings(object, block = block)
    if(block != 0){
      loadingType <- object$info$blockLoadings
      if(!is.null(bn <- names(object$blockLoadings)))
        if(is.numeric(block))
          loadingType <- paste0(loadingType, ", ", bn[block])
      else
        loadingType <- paste0(loadingType, ", ", block)
    }
    else
      loadingType <- object$info$loadings
    if (is.null(L))
      stop("`", deparse(substitute(object)), "' has no loadings.")
  }
  evar <- attr(L,'explvar')[comps]
  L <- L[,comps, drop = FALSE]
  if(is.null(evar))
    varlab <- colnames(L)
  else
    varlab <- paste(colnames(L), " (", format(evar, digits = 2, trim = TRUE),
                    " %)", sep = "")
  if(missing(main))
    main <- ifelse(loadingType != "Not used", loadingType, "")
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
           lwd = lwd, pch = pch, cex = cex, col = col, main = main, ...)
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
      pairs(L, labels = varlab, panel = panel, cex = cex, main = main, ...)
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
            xaxt = xaxt, xlim = xlim, main = main, ...)
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

#' @export
#' @rdname multiblock_plots
loadingweightplot <- function(object, main = "Loading weights", ...){
  mf <- match.call(expand.dots = FALSE)
  object$loadings <- object$loading.weights
  loadingplot(object, main = main, ...)
}


#' @export
#' @rdname multiblock_plots
biplot.multiblock <- function(x, block = 0, comps = 1:2, which = c("x", "y", "scores", "loadings"),
                              var.axes = FALSE, xlabs, ylabs, main, ...)
{
  if (length(comps) != 2) stop("Exactly 2 components must be selected.")
  which <- match.arg(which)
  switch(which,
         x = {
           objects <- scores(x, block = block)
           vars <- loadings(x, block = block)
           title <- "X scores and X loadings"
         },
         y = {
           objects <- x$Yscores
           vars <- x$Yloadings
           title <- "Y scores and Y loadings"
         },
         scores = {
           objects <- scores(x, block = block)
           vars <- x$Yscores
           title <- "X scores and Y scores"
         },
         loadings = {
           objects <- loadings(x, block = block)
           vars <- x$Yloadings
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