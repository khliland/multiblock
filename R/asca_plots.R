#' @title ASCA Result Methods
#' @name asca_plots
#' @aliases asca_plots scoreplot.asca loadingplot.asca
#' 
#' @description Various plotting procedures for \code{\link{asca}} objects. 
#' 
#' @details Usage of the functions are shown using generics in the examples in \code{\link{asca}}.
#' Plot routines are available as
#' \code{scoreplot.asca} and \code{loadingplot.asca}.
#' 
#' @param object \code{asca} object.
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
#' @param ... additional arguments to underlying methods.
#' 
#' @return The plotting routines have no return.
#'
#' @references 
#' * Smilde, A., Jansen, J., Hoefsloot, H., Lamers,R., Van Der Greef, J., and Timmerman, M.(2005). ANOVA-Simultaneous Component Analysis (ASCA): A new tool for analyzing designed metabolomics data. Bioinformatics, 21(13), 3043–3048.
#' * Liland, K.H., Smilde, A., Marini, F., and Næs,T. (2018). Confidence ellipsoids for ASCA models based on multivariate regression theory. Journal of Chemometrics, 32(e2990), 1–13.
#' * Martin, M. and Govaerts, B. (2020). LiMM-PCA: Combining ASCA+ and linear mixed models to analyse high-dimensional designed data. Journal of Chemometrics, 34(6), e3232.
#' 
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results are found in \code{\link{asca_results}}.
#'
#' @export
loadingplot.asca <- function(object, factor = 1, comps = 1:2, ...){
  plot(loadings(object=object, factor=factor), comps=comps, ...)
}

#' @rdname asca_plots
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
  evar <- attr(scors, 'explvar')
  if(missing(xlab))
    xlab <- paste0("Comp ", comps[1], " (",format(evar[comps[1]], digits = 2, trim = TRUE), " %)")
  if(missing(ylab))
    if(length(comps)>1)
      ylab <- paste0("Comp ", comps[2], " (",format(evar[comps[2]], digits = 2, trim = TRUE), " %)")
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
