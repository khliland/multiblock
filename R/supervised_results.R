#' @name supervised_results
#' @title Result functions for Supervised Multiblock Methods
#' @aliases print.mbpls summary.mbpls blockLoadings.mbpls blockScores.mbpls
#' @examples 
#' data(potato)
#' mb <- mbpls(potato[c('Chemical','Compression')], potato[['Sensory']], ncomp = 5)
#' print(mb)
#' 
#' Tb1 <- blockScores(mb)
#' plot(Tb1)
#' 
#' @export
print.mbpls <- function(x, ...){
  cat("Multiblock Partial Least Squares")
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  invisible(x)
}

#' @rdname supervised_results
#' @export
summary.mbpls <- function(object, ...){
  warning('Not implemented yet!')
}

#' @rdname supervised_results
#' @export
blockLoadings <- function (object, ...) {
  UseMethod("blockLoadings", object)
}
#' @rdname supervised_results
#' @export
blockLoadings.mbpls <- function(object, block=1, ...){
  bl <- object$blockLoadings[[block]]
  class(bl) <- "loadings"
  return(bl)
}

#' @rdname supervised_results
#' @export
blockScores <- function (object, ...) {
  UseMethod("blockScores", object)
}
#' @rdname supervised_results
#' @export
blockScores.mbpls <- function(object, block=1, ...){
  sc <- object$blockScores[[block]]
  class(sc) <- "scores"
  return(sc)
}

#' @rdname supervised_results
#' @export
print.mbrda <- function(x, ...){
  cat("Multiblock Redundancy Analysis")
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  invisible(x)
}
