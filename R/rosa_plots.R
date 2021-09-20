#' @name rosa_plots
#' @title Plotting functions for ROSA models
#' 
#' @description Various plotting procedures for \code{\link{rosa}} objects. 
#' 
#' @details Usage of the functions are shown using generics in the examples below. \code{image.rosa}
#' makes an image plot of each candidate score's correlation to the winner or the block-wise
#' response residual. These plots can be used to find alternative block selection for tweaking
#' the ROSA model. \code{barplot.rosa} makes barplot of block and component explained variances.
#' \code{loadingweightsplot} is an adaptation of \code{pls::loadingplot} to plot loading weights.
#'
#' @aliases image.rosa barplot.rosa
#' @param x A \code{rosa} object
#' @param height A \code{rosa} object.
#' @param ncomp Integer to control the number of components to plot (if fewer than the fitted number of components).
#' @param type An optional \code{character} for selecting the plot type. For \code{image.rosa} the options are: "correlation" (default), "residual" or "order". For \code{barplot.rosa} the options indicate: explained variance should be based on training data ("train") or cross-validation ("CV").
#' @param col Colours used for the image and bar plot, defaulting to mcolors(128).
#' @param legend Logical indicating if a legend should be included (default = TRUE) for \code{image.rosa}.
#' @param mar Figure margins, default = c(5,6,4,7) for \code{image.rosa}.
#' @param las Axis text direction, default = 1 for \code{image.rosa}.
#' @param ... Additional parameters passed to \code{loadingplot}, \code{image}, \code{axis}, \code{color.legend}, or \code{barplot}.
#'
#' @return No return.
#' @references Liland, K.H., Næs, T., and Indahl, U.G. (2016). ROSA - a fast extension of partial least squares regression for multiblock data analysis. Journal of Chemometrics, 30, 651–662, doi:10.1002/cem.2824.
#' @importFrom graphics image barplot
#' @importFrom plotrix color.legend
#'
#' @examples
#' data(potato)
#' mod <- rosa(Sensory[,1] ~ ., data = potato, ncomp = 5)
#' image(mod)
#' barplot(mod)
#' loadingweightplot(mod)
#' 
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results in \code{\link{rosa_results}}.
#' @export
image.rosa <- function(x, type = c("correlation","residual","order"), ncomp = x$ncomp,
                       col = mcolors(128), legend = TRUE, mar = c(5,6,4,7), las = 1, ...){
  mf <- match.call(expand.dots = FALSE)
  if(type[1] == "correlation"){
    im <- t(x$candidate.correlation)[1:ncomp,,drop=FALSE]
    if(is.na(match("zlim", names(mf))))
      zlim <- c(-1,1)
    leg <- seq(-1,1,length.out = length(col))
  } else {
    if(type[1] == "residual"){
      im <- -t(x$candidate.RMSE)[1:ncomp,,drop=FALSE]
      if(is.na(match("zlim", names(mf))))
        zlim <- c(min(im),max(im))
      leg <- seq(min(-im),max(-im),length.out = length(col))
    } else {
      if(type[1] == "order"){
        im <- t(x$candidate.RMSE)[1:ncomp,,drop=FALSE]
        for(i in 1:nrow(im)){
          im[i,] <- im[i,]-min(im[i,])
          im[i,] <- 1-im[i,]/max(im[i,])
          if(is.na(match("zlim", names(mf))))
            zlim <- c(min(im),max(im))
        }
        leg <- seq(min(im),max(im),length.out = length(col))
      } else {
        stop("Unsupported plot type")
      }
    }
  }
  legw <- whichMins(legp <- pretty(leg),leg)
  legs <- character(length(col)); legs[legw[!is.na(legw)]] <- legp[!is.na(legw)]

  nresp <- ncol(im)
  if(is.na(match("main", names(mf))))
    main <- ifelse(type[1] == "correlation", "Candidate score correlations", ifelse(type[1] == "residual","Candidate component RMSE","Candidate component residual order"))

  pars <- par(mar = mar, las = las)
  on.exit(par(pars))
  image(im, axes = FALSE, col = col, zlim = zlim, main = main, ...)
  axis(1, at = (0:(ncomp-1))/(ncomp-1), 1:ncomp, ...)
  axis(2, at = (0:(nresp-1))/(nresp-1), colnames(im), ...)
  box()
  for(i in 1:ncomp){
    w <- x$order[[i]]
    points(rep((i-1)/(ncomp-1), length(w)), (w-1)/(nresp-1), pch = 16, col = "white")
    points(rep((i-1)/(ncomp-1), length(w)), (w-1)/(nresp-1))
  }
  if(type[1] == "residual") col <- rev(col)
  color.legend(1.14,0,1.19,1, rect.col = col, gradient = TRUE, legend = legs, ...)
}


#' @export
#' @rdname rosa_plots
barplot.rosa <- function(height, type = c("train","CV"), ncomp = height$ncomp, col = mcolors(ncomp),  ...){
  nums <- attr(numsb <- blockexpl(height, type = type[1], ncomp = ncomp), "compwise")
  nums <- nums[nrow(nums),,drop=FALSE]
  maxCount <- max(height$count)
  nam <- colnames(numsb)[-ncol(numsb)]
  mat <- matrix(0, ncomp, length(nam))
  counts <- rep(1, length(nam))
  corder <- attr(numsb, 'index')
  for(i in 1:ncomp){
    ids <- unlist(lapply(corder, function(j)identical(j,height$order[[i]])))
    mat[i, ids] <- nums[i]
    counts[ids] <- counts[ids] + 1
  }
  colnames(mat) <- nam
  barplot(mat, col = col, ...)
}

