#' @name unsupervised_results
#' @title Result functions for Unsupervised Multiblock Methods
#' @aliases print.sca
#' @examples 
#' data(wine)
#' sc <- sca(wine[c('Smell at rest', 'View', 'Smell after shaking')], ncomp = 4)
#' print(sc)
#' 
#' @export
print.sca <- function(x, ...){
  cat("Simultaneous Component Analysis")
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  invisible(x)
}

#' @rdname unsupervised_results
#' @export
blockScores.sca <- function(object, block=1, ...){
  sc <- object$scores[[block]]
  class(sc) <- "scores"
  return(sc)
}

#' @rdname unsupervised_results
#' @export
blockLoadings.sca <- function(object, block=1, ...){
  bl <- object$loadings[[block]]
  class(bl) <- "loadings"
  return(bl)
}

# TODO:
# pcagca: print, summary, blockLoadings, blockScores