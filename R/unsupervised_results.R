#' @name unsupervised_results
#' @title Result functions for Unsupervised Multiblock Methods
#' @aliases print.sca
#' @examples 
#' data(wine)
#' sc <- sca(wine[c('Smell at rest', 'View', 'Smell after shaking')], ncomp = 4)
#' print(sc)
#' 
#' @export
scores.multiblock <- function(object, block = 0, ...){
  if(block==0 && is.null(object$scores)){
    warning('No global/consensus scores. Returning block 1 scores.')
    block <- 1
  }
  if(block>0 && is.null(object$blockScores)){
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
  if(block>0 && is.null(object$blockLoadings)){
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
  cat(object$info$method, "\n")
  cat("\nCall:\n", deparse(object$call), "\n", sep = "")
  invisible(object)
}

# TODO:
# gca: gca.rgcca?
# pcagca: summary
# mcoa: include or not?
# jive: post-process with PCA?
# All methods: Check that dimnames are given to scores and loadings (- gca.rgcca, mcoa, jive)

# #' @export
# print.sca <- function(x, ...){
#   cat("Simultaneous Component Analysis")
#   cat("\nCall:\n", deparse(x$call), "\n", sep = "")
#   invisible(x)
# }

# #' @rdname unsupervised_results
# #' @export
# blockScores.sca <- function(object, block=1, ...){
#   sc <- object$scores[[block]]
#   class(sc) <- "scores"
#   return(sc)
# }

# #' @rdname unsupervised_results
# #' @export
# blockLoadings.sca <- function(object, block=1, ...){
#   bl <- object$loadings[[block]]
#   class(bl) <- "loadings"
#   return(bl)
# }
