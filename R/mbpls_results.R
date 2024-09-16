#' @title Predict Method for MBPLS
#'
#' @description Prediction for the mbpls (MBPLS) model. New responses or scores are
#' predicted using a fitted model and a data.frame or list containing matrices of observations.
#'
#' @details When \code{type} is \code{"response"} (default), predicted response values
#' are returned.  If \code{comps} is missing (or is \code{NULL}), predictions
#' for \code{length(ncomp)} models with \code{ncomp[1]} components,
#' \code{ncomp[2]} components, etc., are returned.  Otherwise, predictions for
#' a single model with the exact components in \code{comps} are returned.
#' (Note that in both cases, the intercept is always included in the
#' predictions.  It can be removed by subtracting the \code{Ymeans} component
#' of the fitted model.)
#'
#' When \code{type} is \code{"scores"}, predicted score values are returned for
#' the components given in \code{comps}.  If \code{comps} is missing or
#' \code{NULL}, \code{ncomps} is used instead.
#'
#' @param object an \code{mvr} object.  The fitted model
#' @param newdata a data frame.  The new data.  If missing, the training data
#' is used.
#' @param ncomp,comps vector of positive integers.  The components to use in
#' the prediction.  See below.
#' @param type character.  Whether to predict scores or response values
#' @param na.action function determining what should be done with missing
#' values in \code{newdata}.  The default is to predict \code{NA}.  See
#' \code{\link{na.omit}} for alternatives.
#' @param \dots further arguments.  Currently not used
#' @return When \code{type} is \code{"response"}, a three dimensional array of
#' predicted response values is returned.  The dimensions correspond to the
#' observations, the response variables and the model sizes, respectively.
#'
#' When \code{type} is \code{"scores"}, a score matrix is returned.
#' @note A warning message like \samp{'newdata' had 10 rows but variable(s)
#' found have 106 rows} means that not all variables were found in the
#' \code{newdata} data frame.  This (usually) happens if the formula contains
#' terms like \code{yarn$NIR}.  Do not use such terms; use the \code{data}
#' argument instead.  See \code{\link[pls]{mvr}} for details.
#' @author Kristian Hovde Liland
#' @seealso \code{\link{mbpls}}
#' @examples
#' data(potato)
#' mb <- mbpls(Sensory ~ Chemical+Compression, data=potato, ncomp = 5, subset = 1:26 <= 18)
#' testdata <- subset(potato, 1:26 > 18)
#' 
#' # Predict response
#' yhat <- predict(mb, newdata = testdata)
#' 
#' # Predict scores and plot
#' scores <- predict(mb, newdata = testdata, type = "scores")
#' scoreplot(mb)
#' points(scores[,1], scores[,2], col="red")
#' legend("topright", legend = c("training", "test"), col=1:2, pch = 1)
#' @export
predict.mbpls <- function(object, newdata, ncomp = 1:object$ncomp, comps, 
                          type = c("response", "scores"), 
                          na.action = na.pass, ...){
  class(object) <- c('multiblock','mvr')
  if(!missing(newdata)){
    newdata <- newdata[names(object$blockScores)]
    nblock <- length(newdata)
    Xo <- lapply(newdata, function(x)if(is.factor(x)){return(dummycode(x))}else{return(x)})
    scale <- ifelse(is.null(object$scale), FALSE, object$scale)
    Xo <- lapply(lapply(Xo, as.matrix), function(x)scale(x, scale=scale))
    if(is.numeric(object$blockScale)){ # User supplied scale
      X <- lapply(1:nblock, function(i) Xo[[i]]/norm(Xo[[i]],"F")^2*object$blockScale[i])
    } else {
      if(is.character(object$blockScale)){
        if(object$blockScale[1] == "sqrtnvar")
          X <- lapply(Xo, function(x) x/sqrt(ncol(x)))
        else if(object$blockScale[1] == "ssq")
          X <- lapply(1:nblock, function(i) Xo[[i]]/norm(Xo[[i]],"F")^2)
        else if(object$blockScale[1] == "none")
          X <- Xo
        else
          stop("Unknown format for 'object$blockScale'")
      } else {
        stop("Unknown format for 'object$blockScale'")
      }
    }
    newdata <- do.call(cbind, X)
  }
  return(predict(object = object, newdata = newdata, ncomp = ncomp, comps = comps, 
                 type = type, na.action = na.action, ...))
}

### mvrVal.R: Functions for calculating validation statistics, such
### as MSEP, RMSEP and R2, for mvr objects.

## Calculate the validation statistics needed for (R)MSEP and R^2.
## Note that it accepts any values for `estimate', but only calculates
## statistics for "train", "test" and "CV".

#' @name mvrVal
#' @title MSEP, RMSEP and R2 of the MB-PLS model
#'
#' @description Functions to estimate the mean squared error of prediction (MSEP), root mean
#' squared error of prediction (RMSEP) and \eqn{R^2} (A.K.A. coefficient of
#' multiple determination) for a fitted MB-PLS models.  Test-set,
#' cross-validation and calibration-set estimates are implemented.
#'
#' @details \code{RMSEP} simply calls \code{MSEP} and takes the square root of the
#' estimates.  It therefore accepts the same arguments as \code{MSEP}.
#'
#' Several estimators can be used.  \code{"train"} is the training or
#' calibration data estimate, also called (R)MSEC.  For \code{R2}, this is the
#' unadjusted \eqn{R^2}.  It is overoptimistic and should not be used for
#' assessing models.  \code{"CV"} is the cross-validation estimate, and
#' \code{"adjCV"} (for \code{RMSEP} and \code{MSEP}) is the bias-corrected
#' cross-validation estimate.  They can only be calculated if the model has
#' been cross-validated.  Finally, \code{"test"} is the test set estimate,
#' using \code{newdata} as test set.
#'
#' Which estimators to use is decided as follows (see below for
#' \code{pls:mvrValstats}).  If \code{estimate} is not specified, the test set
#' estimate is returned if \code{newdata} is specified, otherwise the CV and
#' adjusted CV (for \code{RMSEP} and \code{MSEP}) estimates if the model has
#' been cross-validated, otherwise the training data estimate.  If
#' \code{estimate} is \code{"all"}, all possible estimates are calculated.
#' Otherwise, the specified estimates are calculated.
#'
#' Several model sizes can also be specified.  If \code{comps} is missing (or
#' is \code{NULL}), \code{length(ncomp)} models are used, with \code{ncomp[1]}
#' components, \ldots{}, \code{ncomp[length(ncomp)]} components.  Otherwise, a
#' single model with the components \code{comps[1]}, \ldots{},
#' \code{comps[length(comps)]} is used.  If \code{intercept} is \code{TRUE}, a
#' model with zero components is also used (in addition to the above).
#'
#' The \eqn{R^2} values returned by \code{"R2"} are calculated as \eqn{1 -
#' SSE/SST}, where \eqn{SST} is the (corrected) total sum of squares of the
#' response, and \eqn{SSE} is the sum of squared errors for either the fitted
#' values (i.e., the residual sum of squares), test set predictions or
#' cross-validated predictions (i.e., the \eqn{PRESS}).  For \code{estimate =
#' "train"}, this is equivalent to the squared correlation between the fitted
#' values and the response.  For \code{estimate = "train"}, the estimate is
#' often called the prediction \eqn{R^2}.
#'
#' \code{mvrValstats} is a utility function that calculates the statistics
#' needed by \code{MSEP} and \code{R2}.  It is not intended to be used
#' interactively.  It accepts the same arguments as \code{MSEP} and \code{R2}.
#' However, the \code{estimate} argument must be specified explicitly: no
#' partial matching and no automatic choice is made.  The function simply
#' calculates the types of estimates it knows, and leaves the other untouched.
#'
#' @aliases MSEP.mbpls RMSEP.mbpls R2.mbpls mvrValstats.mbpls
#' @param object an \code{mvr} object
#' @param estimate a character vector.  Which estimators to use.  Should be a
#' subset of \code{c("all", "train", "CV", "adjCV", "test")}.  \code{"adjCV"}
#' is only available for (R)MSEP.  See below for how the estimators are chosen.
#' @param newdata a data frame with test set data.
#' @param ncomp,comps a vector of positive integers.  The components or number
#' of components to use.  See below.
#' @param intercept logical.  Whether estimates for a model with zero
#' components should be returned as well.
#' @param se logical.  Whether estimated standard errors of the estimates
#' should be calculated.  Not implemented yet.
#' @param \dots further arguments sent to underlying functions or (for
#' \code{RMSEP}) to \code{MSEP}
#' @section Value: \code{mvrValstats} returns a list with components \describe{
#' \item{SSE}{three-dimensional array of SSE values.  The first dimension is
#' the different estimators, the second is the response variables and the third
#' is the models.} \item{SST}{matrix of SST values.  The first dimension is the
#' different estimators and the second is the response variables.}
#' \item{nobj}{a numeric vector giving the number of objects used for each
#' estimator.} \item{comps}{the components specified, with \code{0} prepended
#' if \code{intercept} is \code{TRUE}.} \item{cumulative}{\code{TRUE} if
#' \code{comps} was \code{NULL} or not specified.} }
#'
#' The other functions return an object of class \code{"mvrVal"}, with
#' components \describe{ \item{val}{three-dimensional array of estimates.  The
#' first dimension is the different estimators, the second is the response
#' variables and the third is the models.} \item{type}{\code{"MSEP"},
#' \code{"RMSEP"} or \code{"R2"}.} \item{comps}{the components specified, with
#' \code{0} prepended if \code{intercept} is \code{TRUE}.}
#' \item{cumulative}{\code{TRUE} if \code{comps} was \code{NULL} or not
#' specified.} \item{call}{the function call} }
#' @author Kristian Hovde Liland
#' @seealso \code{\link{mbpls}}
#' @references Mevik, B.-H., Cederkvist, H. R. (2004) Mean Squared Error of
#' Prediction (MSEP) Estimates for Principal Component Regression (PCR) and
#' Partial Least Squares Regression (PLSR).  \emph{Journal of Chemometrics},
#' \bold{18}(9), 422--429.
#' @keywords regression multivariate
#' @examples
#'
#' data(oliveoil, package = "pls")
#' mod <- pls::plsr(sensory ~ chemical, ncomp = 4, data = oliveoil, validation = "LOO")
#' RMSEP(mod)
#' \dontrun{plot(R2(mod))}
#'
#' @export
R2.mbpls <- function(object, estimate, newdata, ncomp = 1:object$ncomp, comps,
                     intercept = TRUE, se = FALSE, ...)
{
  cumulative <- missing(comps) || is.null(comps)
  class(object) <- c('multiblock','mvr')
  if(!missing(newdata)){
    newdata <- newdata[names(object$blockScores)]
    nblock <- length(newdata)
    Xo <- lapply(newdata, function(x)if(is.factor(x)){return(dummycode(x))}else{return(x)})
    scale <- ifelse(is.null(object$scale), FALSE, object$scale)
    Xo <- lapply(lapply(Xo, as.matrix), function(x)scale(x, scale=scale))
    if(is.numeric(object$blockScale)){ # User supplied scale
      X <- lapply(1:nblock, function(i) Xo[[i]]/norm(Xo[[i]],"F")^2*object$blockScale[i])
    } else {
      if(is.character(object$blockScale)){
        if(object$blockScale[1] == "sqrtnvar")
          X <- lapply(Xo, function(x) x/sqrt(ncol(x)))
        else if(object$blockScale[1] == "ssq")
          X <- lapply(1:nblock, function(i) Xo[[i]]/norm(Xo[[i]],"F")^2)
        else if(object$blockScale[1] == "none")
          X <- Xo
        else
          stop("Unknown format for 'object$blockScale'")
      } else {
        stop("Unknown format for 'object$blockScale'")
      }
    }
    newdata$X <- I(do.call(cbind, X))
    browser()
    respnames <- dimnames(fitted(object))[[2]]
  }
  return(R2(object = object, estimate = estimate, newdata = newdata, 
            ncomp = ncomp, comps=comps, intercept = intercept, se = se, ...))
}


## MSEP: Return MSEP
#' @rdname mvrVal
#' @export
MSEP.mbpls <- function(object, estimate, newdata, ncomp = 1:object$ncomp, comps,
                       intercept = TRUE, se = FALSE, ...)
{
  class(object) <- c('multiblock','mvr')
  if(!missing(newdata)){
    newdata <- newdata[names(object$blockScores)]
    nblock <- length(newdata)
    Xo <- lapply(newdata, function(x)if(is.factor(x)){return(dummycode(x))}else{return(x)})
    scale <- ifelse(is.null(object$scale), FALSE, object$scale)
    Xo <- lapply(lapply(Xo, as.matrix), function(x)scale(x, scale=scale))
    if(is.numeric(object$blockScale)){ # User supplied scale
      X <- lapply(1:nblock, function(i) Xo[[i]]/norm(Xo[[i]],"F")^2*object$blockScale[i])
    } else {
      if(is.character(object$blockScale)){
        if(object$blockScale[1] == "sqrtnvar")
          X <- lapply(Xo, function(x) x/sqrt(ncol(x)))
        else if(object$blockScale[1] == "ssq")
          X <- lapply(1:nblock, function(i) Xo[[i]]/norm(Xo[[i]],"F")^2)
        else if(object$blockScale[1] == "none")
          X <- Xo
        else
          stop("Unknown format for 'object$blockScale'")
      } else {
        stop("Unknown format for 'object$blockScale'")
      }
    }
    newdata <- do.call(cbind, X)
  }
  return(MSEP(object = object, estimate = estimate, newdata = newdata, 
              ncomp = ncomp, comps = comps, intercept = intercept, se = se, ...))
}

# RMSEP: A wrapper around MSEP to calculate RMSEPs
#' @rdname mvrVal
#' @export
RMSEP.mbpls <- function(object, ...) {
  class(object) <- c('multiblock','mvr')
  if(!missing(newdata)){
    newdata <- newdata[names(object$blockScores)]
    nblock <- length(newdata)
    Xo <- lapply(newdata, function(x)if(is.factor(x)){return(dummycode(x))}else{return(x)})
    scale <- ifelse(is.null(object$scale), FALSE, object$scale)
    Xo <- lapply(lapply(Xo, as.matrix), function(x)scale(x, scale=scale))
    if(is.numeric(object$blockScale)){ # User supplied scale
      X <- lapply(1:nblock, function(i) Xo[[i]]/norm(Xo[[i]],"F")^2*object$blockScale[i])
    } else {
      if(is.character(object$blockScale)){
        if(object$blockScale[1] == "sqrtnvar")
          X <- lapply(Xo, function(x) x/sqrt(ncol(x)))
        else if(object$blockScale[1] == "ssq")
          X <- lapply(1:nblock, function(i) Xo[[i]]/norm(Xo[[i]],"F")^2)
        else if(object$blockScale[1] == "none")
          X <- Xo
        else
          stop("Unknown format for 'object$blockScale'")
      } else {
        stop("Unknown format for 'object$blockScale'")
      }
    }
    newdata <- do.call(cbind, X)
  }
  return(RMSEP(object, ...))
}
