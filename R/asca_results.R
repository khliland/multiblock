#' @title ASCA Result Methods
#' @name asca_object
#' @aliases asca_object print.asca summary.asca projections projections.asca print.summary.asca loadings.asca scores.asca scoreplot.asca
#' @param object \code{asca} object.
#' @param x \code{asca} object.
#' @param factor \code{integer/character} for selecting a model factor.
#' @param comps \code{integer} vector of selected components.
#' @param pch.scores \code{integer} plotting symbol.
#' @param pch.projections \code{integer} plotting symbol.
#' @param gr.col \code{integer} vector of colours for groups.
#' @param ellipsoids \code{character} "confidence" or "data" ellipsoids for balanced fixed effect models.
#' @param xlim \code{numeric} x limits.
#' @param ylim \code{numeric} y limits. 
#' @param xlab \code{character} x label.
#' @param ylab \code{character} y label. 
#' @param legendpos \code{character} position of legend.
#' @param digits \code{integer} number of digits for printing.
#' @param ... addtitional arguments to underlying methods.
#'
#' @export
print.asca <- function(x, ...){
  cat("Anova Simultaneous Component Analysis fitted using", x$fit.type)
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  invisible(x)
}

#' @rdname asca_object
#' @export
summary.asca <- function(object, ...){
  dat <- data.frame(ss=object$ssq, expl=object$explvar)
  dat <- dat[-nrow(dat),,drop=FALSE]
  x <- list(dat=dat, fit.type=object$fit.type)
  class(x) <- c('summary.asca')
  x
}

#' @rdname asca_object
#' @export
print.summary.asca <- function(x, digits=2, ...){
  cat("Anova Simultaneous Component Analysis fitted using", x$fit.type, "\n")
  print(round(x$dat, digits))
  invisible(x$dat)
}

#' @rdname asca_object
#' @export
loadings.asca <- function(object, factor = 1, ...){
  loads <- object$loadings[[factor]]
  class(loads) <- "loadings"
  return(loads)
}
#' @rdname asca_object
#' @export
loadingplot.asca <- function(object, factor = 1, comps = 1:2, ...){
  plot(loadings(object=object, factor=factor), comps=comps, ...)
}

#' @rdname asca_object
#' @export
scores.asca <- function(object, factor = 1, ...){
  scors <- object$scores[[factor]]
  class(scors) <- "scores"
  return(scors)
}

#' @rdname asca_object
#' @export
projections <- function (object, ...) {
  UseMethod("projections", object)
}

#' @rdname asca_object
#' @export
projections.asca <- function(object, factor = 1, ...){
  projs <- object$projected[[factor]]
  class(projs) <- "projs"
  return(projs)
}

#' @rdname asca_object
#' @export
scoreplot.asca <- function(object, factor = 1, comps = 1:2, pch.scores = 19, pch.projections = 1, 
                           gr.col = 1:nlevels(object$effects[[factor]]), ellipsoids,
                           xlim,ylim, xlab,ylab, legendpos, ...){
  # Number of levels in current factor
  nlev  <- nlevels(object$effects[[factor]])
  nobj  <- nrow(object$Y)
  # Remove redundant levels
  comps <- comps[comps <= nlev-1]
  
  scors <- scores(object=object, factor=factor)
  projs <- projections(object=object, factor=factor) + scors
  if(missing(xlim))
    xlim <- c(min(min(scors[,comps[1]]), min(projs[,comps[1]])),
              max(max(scors[,comps[1]]), max(projs[,comps[1]])))
  if(missing(ylim))
    if(length(comps)>1)
      ylim <- c(min(min(scors[,comps[2]]), min(projs[,comps[2]])),
                max(max(scors[,comps[2]]), max(projs[,comps[2]])))
  else
    ylim <- c(0.5, nlev+0.5)
  if(missing(xlab))
    xlab <- paste("Comp", comps[1])
  if(missing(ylab))
    if(length(comps)>1)
      ylab <- paste("Comp", comps[2])
  else
    ylab <- 'Level'
  
  if(length(comps)>1){ # Scatter plot
    scoreplot(scors, comps=comps, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, pch=pch.scores, ...)
    for(i in 1:nlev){
      lev <- levels(object$effects[[factor]])[i]
      points(scors[object$effects[[factor]] == lev, comps], pch=pch.scores, col=gr.col[i])
      points(projs[object$effects[[factor]] == lev, comps], pch=pch.projections, col=gr.col[i])
    }
    if(!missing(legendpos))
      legend(legendpos, legend = levels(object$effects[[factor]]), col=gr.col, pch=pch.scores)
    
    if(!missing(ellipsoids)){
      if(ellipsoids == "data"){
        dataEllipse(projs[,comps], groups = object$effects[[factor]], levels=c(0.4,0.68,0.95), add=TRUE, plot.points=FALSE, col=gr.col, lwd=1, group.labels="", center.pch=FALSE, lty=3)
      }
      if(ellipsoids == "confidence" || ellipsoids == "conf"){
        # Covariance matrix
        sigma <- crossprod(object$residuals)/nobj
        L <- object$loadings[[factor]][,comps]
        # Transformed covariance matrix
        LSL <- crossprod(L,sigma) %*% L * nlev / (nobj-nlev)
        if(!isSymmetric(LSL))
          LSL <- (LSL + t(LSL))/2 # Force symmetry
        # Scaling by confidence
        c40 <- sqrt((nobj-nlev)*2 / (nobj-nlev-2+1) * qf(0.40, 2, nobj-nlev-2+1))
        c68 <- sqrt((nobj-nlev)*2 / (nobj-nlev-2+1) * qf(0.68, 2, nobj-nlev-2+1))
        c95 <- sqrt((nobj-nlev)*2 / (nobj-nlev-2+1) * qf(0.95, 2, nobj-nlev-2+1))
        for(i in 1:nlev){
          lev <- levels(object$effects[[factor]])[i]
          ellipse(colMeans(scors[object$effects[[factor]]==lev,comps]), LSL, c40, lwd=1, col=gr.col[i])
          ellipse(colMeans(scors[object$effects[[factor]]==lev,comps]), LSL, c68, lwd=1, col=gr.col[i])
          ellipse(colMeans(scors[object$effects[[factor]]==lev,comps]), LSL, c95, lwd=1, col=gr.col[i])
        }
      }
    }
  } else { # Line plot
    plot(scors[,comps], as.numeric(object$effects[[factor]]), xlim=xlim, 
         ylim=ylim, xlab=xlab, ylab=ylab, axes = FALSE)
    axis(1)
    axis(2, at=1:nlev, labels = levels(object$effects[[factor]]))
    box()
    for(i in 1:nlev){
      lev <- levels(object$effects[[factor]])[i]
      points(scors[object$effects[[factor]] == lev, comps], rep(i,sum(as.numeric(object$effects[[factor]]) == i)), pch=pch.scores, col=gr.col[i])
      points(projs[object$effects[[factor]] == lev, comps], rep(i,sum(as.numeric(object$effects[[factor]]) == i)), pch=pch.projections, col=gr.col[i])
    }
    if(!missing(ellipsoids)){
      if(ellipsoids == "confidence" || ellipsoids == "conf"){
        sigma <- crossprod(object$residuals)/nobj
        L <- object$loadings[[factor]][,comps]
        # Transformed covariance matrix
        LSL <- sqrt(crossprod(L,sigma) %*% L * nlev / (nobj-nlev))
        # Scaling by confidence
        c40 <- sqrt((nobj-nlev)*1 / (nobj-nlev-1+1) * qf(0.40, 1, nobj-nlev-1+1))
        c68 <- sqrt((nobj-nlev)*1 / (nobj-nlev-1+1) * qf(0.68, 1, nobj-nlev-1+1))
        c95 <- sqrt((nobj-nlev)*1 / (nobj-nlev-1+1) * qf(0.95, 1, nobj-nlev-1+1))
        for(i in 1:nlev){
          lev <- levels(object$effects[[factor]])[i]
          lines(mean(scors[object$effects[[factor]] == lev,comps])*c(1,1)+c(LSL)*c40, i+c(-0.2,0.2), col=gr.col[i])
          lines(mean(scors[object$effects[[factor]] == lev,comps])*c(1,1)+c(LSL)*c68, i+c(-0.2,0.2), col=gr.col[i])
          lines(mean(scors[object$effects[[factor]] == lev,comps])*c(1,1)+c(LSL)*c95, i+c(-0.2,0.2), col=gr.col[i])
          lines(mean(scors[object$effects[[factor]] == lev,comps])*c(1,1)-c(LSL)*c40, i+c(-0.2,0.2), col=gr.col[i])
          lines(mean(scors[object$effects[[factor]] == lev,comps])*c(1,1)-c(LSL)*c68, i+c(-0.2,0.2), col=gr.col[i])
          lines(mean(scors[object$effects[[factor]] == lev,comps])*c(1,1)-c(LSL)*c95, i+c(-0.2,0.2), col=gr.col[i])
        }
      }
    }
  }
}
