#' @title ASCA Result Methods
#' @name asca_results
#' @aliases asca_results print.asca summary.asca projections projections.asca print.summary.asca loadings.asca scores.asca
#' 
#' @description Standard result computation and extraction functions for ASCA (\code{\link{asca}}).
#' 
#' @details Usage of the functions are shown using generics in the examples in \code{\link{asca}}.
#' Explained variances are available (block-wise and global) through \code{blockexpl} and \code{print.rosaexpl}.
#' Object printing and summary are available through: 
#' \code{print.asca} and \code{summary.asca}.
#' Scores and loadings have their own extensions of \code{scores()} and \code{loadings()} through
#' \code{scores.asca} and \code{loadings.asca}. Special to ASCA is that scores are on a
#' factor level basis, while back-projected samples have their own function in \code{projections.asca}.
#' 
#' @param object \code{asca} object.
#' @param x \code{asca} object.
#' @param factor \code{integer/character} for selecting a model factor.
#' @param digits \code{integer} number of digits for printing.
#' @param ... additional arguments to underlying methods.
#'
#' @return Returns depend on method used, e.g. \code{projections.asca} returns projected samples, 
#' \code{scores.asca} return scores, while print and summary methods return the object invisibly.
#' 
#' @references 
#' * Smilde, A., Jansen, J., Hoefsloot, H., Lamers,R., Van Der Greef, J., and Timmerman, M.(2005). ANOVA-Simultaneous Component Analysis (ASCA): A new tool for analyzing designed metabolomics data. Bioinformatics, 21(13), 3043–3048.
#' * Liland, K.H., Smilde, A., Marini, F., and Næs,T. (2018). Confidence ellipsoids for ASCA models based on multivariate regression theory. Journal of Chemometrics, 32(e2990), 1–13.
#' * Martin, M. and Govaerts, B. (2020). LiMM-PCA: Combining ASCA+ and linear mixed models to analyse high-dimensional designed data. Journal of Chemometrics, 34(6), e3232.
#' 
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for plotting are found in \code{\link{asca_plots}}.
#'
#' @export
print.asca <- function(x, ...){
  cat("Anova Simultaneous Component Analysis fitted using", x$fit.type)
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  invisible(x)
}

#' @rdname asca_results
#' @export
summary.asca <- function(object, ...){
  dat <- data.frame(ss=object$ssq, expl=object$explvar)
  dat <- dat[-nrow(dat),,drop=FALSE]
  x <- list(dat=dat, fit.type=object$fit.type)
  class(x) <- c('summary.asca')
  x
}

#' @rdname asca_results
#' @export
print.summary.asca <- function(x, digits=2, ...){
  cat("Anova Simultaneous Component Analysis fitted using", x$fit.type, "\n")
  print(round(x$dat, digits))
  invisible(x$dat)
}

#' @rdname asca_results
#' @export
loadings.asca <- function(object, factor = 1, ...){
  loads <- object$loadings[[factor]]
  class(loads) <- "loadings"
  return(loads)
}

#' @rdname asca_results
#' @export
scores.asca <- function(object, factor = 1, ...){
  scors <- object$scores[[factor]]
  class(scors) <- "scores"
  return(scors)
}

#' @rdname asca_results
#' @export
projections <- function (object, ...) {
  UseMethod("projections", object)
}

#' @rdname asca_results
#' @export
projections.asca <- function(object, factor = 1, ...){
  projs <- object$projected[[factor]]
  class(projs) <- "projs"
  return(projs)
}
