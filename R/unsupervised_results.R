#' @name unsupervised_results
#' @title Result functions for Unsupervised Multiblock Methods
#' @aliases print.sca
#' 
#' @param object \code{multiblock} object.
#' @param x \code{multiblock} object.
#' @param block \code{integer/character} for block selection.
#' @param ... Not implemented.

#' @examples 
#' data(wine)
#' sc <- sca(wine[c('Smell at rest', 'View', 'Smell after shaking')], ncomp = 4)
#' print(sc)
#' plot(loadings(sc, block = 1), labels = "names", scatter = TRUE)
#' 
#' @seealso Overviews of available methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' @export
scores.multiblock <- function(object, block = 0, ...){
  if(block==0 && is.null(object$scores)){
    warning('No global/consensus scores. Returning block 1 scores.')
    block <- 1
  }
  if(block!=0 && is.null(object$blockScores)){
    warning('No block scores. Returning global/consensus scores.')
    block <- 0
  }
  if(block==0){
    s <- object$scores
    attr(s, 'info') <- object$info$scores
  } else {
    s <- object$blockScores[[block]]
    attr(s, 'info') <- object$info$blockScores
  }
  class(s) <- list('scores')
  return(s)
}

#' @rdname unsupervised_results 
#' @export
loadings.multiblock <- function(object, block = 0, ...){
  if(block==0 && is.null(object$loadings)){
    warning('No global/consensus loadings available. Returning block 1 loadings.')
    block <- 1
  }
  if(block!=0 && is.null(object$blockLoadings)){
    warning('No block loadings available. Returning global/consensus loadings.')
    block <- 0
  }
  if(block==0){
    l <- object$loadings
    attr(l, 'info') <- object$info$loadings
  } else {
    l <- object$blockLoadings[[block]]
    attr(l, 'info') <- object$info$blockLoadings
  }
  class(l) <- list('loadings')
  return(l)
}

#' @rdname unsupervised_results 
#' @export
print.multiblock <- function(x, ...){
  cat(x$info$method, "\n")
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  invisible(x)
}

#' @rdname unsupervised_results 
#' @export
summary.multiblock <- function(object, ...){
  cat(object$info$method, "\n\n")
  if(!is.null(object$scores))
    cat("Score type: ", object$info$scores, "\n")
  if(!is.null(object$loadings))
    cat("Loadings type: ", object$info$loadings, "\n")
  if(!is.null(object$blockScores))
    cat("Block scores type: ", object$info$blockScores, "\n")
  if(!is.null(object$blockLoadings))
    cat("Block loadings type: ", object$info$blockLoadings, "\n")
  invisible(object)
}

# #' @rdname unsupervised_results
# #' @export
# info <- function (object, ...) {
#   UseMethod("info", object)
# }
# 
# #' @rdname unsupervised_results
# #' @export
# info.multiblock(object, ...) {
#   summary(object)
# }
