#' @name multiblock_results
#' @title Result Functions for Multiblock Objects
#' @aliases scores.multiblock loadings.multiblock print.multiblock summary.multiblock
#' 
#' @description Standard result computation and extraction functions for \code{multiblock} objects.
#' 
#' @details Usage of the functions are shown using generics in the examples below.
#' Object printing and summary are available through: 
#' \code{print.multiblock} and \code{summary.multiblock}.
#' Scores and loadings have their own extensions of \code{scores()} and \code{loadings()} throught
#' \code{scores.multiblock} and \code{loadings.multiblock}.
#' 
#' @param object \code{multiblock} object.
#' @param x \code{multiblock} object.
#' @param block \code{integer/character} for block selection.
#' @param ... Not implemented.
#' 
#' @return Scores or loadings are returned by \code{scores.multiblock} and \code{loadings.multiblock}, while print and summary methods invisibly returns the object.
#' 
#' @examples 
#' data(wine)
#' sc <- sca(wine[c('Smell at rest', 'View', 'Smell after shaking')], ncomp = 4)
#' print(sc)
#' summary(sc)
#' head(loadings(sc, block = 1))
#' head(scores(sc))
#' 
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for plotting are found in \code{\link{multiblock_plots}}, respectively.
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
  class(s) <- list('scores.multiblock','scores')
  return(s)
}

#' @rdname multiblock_results 
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
  class(l) <- list('loadings.multiblock','loadings')
  return(l)
}

#' @rdname multiblock_results 
#' @export
print.multiblock <- function(x, ...){
  cat(x$info$method, "\n")
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  invisible(x)
}

#' @rdname multiblock_results 
#' @export
summary.multiblock <- function(object, ...){
  cat(object$info$method, "\n")
  cat(paste0(rep("=", nchar(object$info$method)), collapse=""), "\n\n")
  if(!is.null(object$scores)){
    cat("$scores: ", object$info$scores, " (",nrow(object$scores),"x",ncol(object$scores),")", "\n", sep="")
  }
  if(!is.null(object$loadings)){
    cat("$loadings: ", object$info$loadings, " (",nrow(object$loadings),"x",ncol(object$loadings),")", "\n", sep="")
  }
  if(!is.null(object$blockScores)){
    cat("$blockScores: ", object$info$blockScores, ":\n", sep="")
    bs <- lapply(lapply(object$blockScores, dim), function(x)paste0("(",x[1],"x",x[2],")"))
    cat("- ", paste(apply(cbind(names(bs), unlist(bs)), 1, paste, collapse=" "), collapse = ", "), "\n", sep="")
  }
  if(!is.null(object$blockLoadings)){
    cat("$blockLoadings: ", object$info$blockLoadings, ":\n", sep="")
    bl <- lapply(lapply(object$blockLoadings, dim), function(x)paste0("(",x[1],"x",x[2],")"))
    cat("- ", paste(apply(cbind(names(bl), unlist(bl)), 1, paste, collapse=" "), collapse = ", "), "\n", sep="")
  }
  invisible(object)
  # Extend with dimensions and number of blocks
}
