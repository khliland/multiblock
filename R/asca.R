#' @name asca
#' @aliases asca print.asca summary.asca loadings.asca scores.asca loadingplot.asca scoreplot.asca
#' @title Analysis of Variance Simultaneous Component Analysis
#'
#' @param formula Model formula accepting a single response (block) and predictor names separated by + signs.
#' @param data The data set to analyse.
#' @param subset Expression for subsetting the data before modelling.
#' @param weights Optional object weights.
#' @param subset Subset of objects
#' @param na.action How to handle NAs (no action implemented).
#' @param family Error distributions and link function for Generalized Linear Models.
#' @param pca.in Compress response before ASCA (number of components).
#' @param object \code{ASCA} object used when extracting and plotting.
#' @param factor String or integer indicating which factor to select for extracting and plotting.
#'
#' @return
#'
#' @importfrom lme4 lmer
#' @examples
#' dataset   <- data.frame(y=I(matrix(rnorm(24*10),ncol=10)), x=factor(c(rep(2,8),rep(1,8),rep(0,8))), z=factor(rep(c(1,0),12)), w=rnorm(24))
#' colnames(dataset$y) <- paste('Var', 1:10, sep=" ")
#' rownames(dataset) <- paste('Obj', 1:24, sep=" ")
#' mod <- asca(y~x+z, data=dataset)
#' print(mod)
#' loadingplot(mod, scatter=TRUE)
#' scoreplot(mod)
#' 
#' mod.pca <- asca(y~x+z, data=dataset, pca.in=5)
#' mod.glm <- asca(y~x+z+w, data=dataset, family="gaussian")
#' 
#' mod <- asca(y~x+(1|z), data=dataset)
#' mod <- asca(y~x+(1|z)+w, data=dataset, family="gaussian")
#' @export
asca <- function(formula, data, subset, weights, na.action, family, pca.in = FALSE){
  ## Force contrast to sum
  opt <- options(contrasts = c(unordered="contr.sum", ordered="contr.poly"))
  
  ## Get the data matrices
  Y <- data[[formula[[2]]]]
  N <- nrow(Y)
  Y <- Y - rep(colMeans(Y), each=N) # Centre Y
  ssqY <- sum(Y^2)
  if(pca.in != 0){
    if(pca.in == 1)
      stop('pca.in = 1 is not supported (single response)')
    Yudv <- svd(Y)
    Y <- Yudv$u[,1:pca.in,drop=FALSE] * rep(Yudv$d[1:pca.in], each=N)
  }
  residuals <- Y
  
  mf <- match.call(expand.dots = FALSE)
  fit.type <- "'lm' (Linear Model)"
  if(length(grep('|', formula, fixed=TRUE)) == 0){
    # Fixed effect model
    if(missing(family)){
      # LM
      m <- match(c("formula", "data", "weights", "subset", "na.action"), names(mf), 0)
      mf <- mf[c(1, m)]                # Retain only the named arguments
      mf[[1]] <- as.name("lm")
      mf[[3]] <- as.name("dat")
      dat <- data
      dat[[formula[[2]]]] <- Y
      ano   <- eval(mf, envir = environment())
      coefs <- as.matrix(coefficients(ano))
    } else {
      # GLM
      m <- match(c("formula", "data", "weights", "subset", "na.action", "family"), names(mf), 0)
      mf <- mf[c(1, m)]                # Retain only the named arguments
      mf[[1]] <- as.name("glm")
      mf[[3]] <- as.name("dat")
      dat <- data
      for(i in 1:ncol(Y)){
        dat[[formula[[2]]]] <- Y[,i,drop=FALSE]
        ano <- eval(mf, envir = environment())
        if(i == 1)
          coefs <- matrix(0.0, length(coefficients(ano)), ncol(Y))
        coefs[,i] <- coefficients(ano)
      }
      fit.type <- "'glm' (Generalized Linear Model)"
    }
  } else {
    # Mixed model
    if(missing(family)){
      # LM
      m <- match(c("formula", "data", "weights", "subset", "na.action"), names(mf), 0)
      mf <- mf[c(1, m)]                # Retain only the named arguments
      mf[[1]] <- as.name("lmer")
      mf[[3]] <- as.name("dat")
      dat <- data
      for(i in 1:ncol(Y)){
        dat[[formula[[2]]]] <- Y[,i,drop=FALSE]
        ano <- eval(mf, envir = environment())
        if(i == 1)
          coefs <- matrix(0.0, length(colMeans(coefficients(ano)[[1]])), ncol(Y))
        coefs[,i] <- colMeans(coefficients(ano)[[1]])
      }
      fit.type <- "'lmer' (Linear Mixed Model)"
    } else {
      # GLM
      m <- match(c("formula", "data", "weights", "subset", "na.action", "family"), names(mf), 0)
      mf <- mf[c(1, m)]                # Retain only the named arguments
      mf[[1]] <- as.name("glmer")
      mf[[3]] <- as.name("dat")
      dat <- data
      for(i in 1:ncol(Y)){
        dat[[formula[[2]]]] <- Y[,i,drop=FALSE]
        ano <- eval(mf, envir = environment())
        if(i == 1) # colMeans assumes only random intercepts, not slopes
          coefs <- matrix(0.0, length(colMeans(coefficients(ano)[[1]])), ncol(Y))
        coefs[,i] <- colMeans(coefficients(ano)[[1]])
      }
      fit.type <- "'glmer' (Generalized Linear Mixed Model)"
    }
  }
  M      <- model.matrix(ano)
  effs   <- attr(terms(ano), "term.labels")
  assign <- attr(M, "assign")
  modFra <- model.frame(ano)
  
  # Exclude numeric effects and their interactions
  nums   <- names(which(unlist(lapply(modFra, class)) == "numeric"))
  if(length(nums)>0){
    exclude  <- match(nums, rownames(attr(terms(ano), "factors")))
    approved <- which(colSums(attr(terms(ano), "factors")[exclude,,drop=FALSE])==0)
  } else {
    approved <- 1:max(assign)
  }
  LS <- effects <- ssq <- list()
  for(i in approved){
    LS[[effs[i]]] <- M[, assign==i, drop=FALSE] %*% coefs[assign==i,]
    residuals <- residuals - LS[[effs[i]]]
    effects[[effs[i]]] <- modFra[[effs[i]]]
    ssq[[effs[i]]] <- sum(LS[[effs[i]]]^2)
  }
  ssq$res <- sum(residuals^2)

  # SCAs
  scores <- loadings <- projected <- singulars <- list()
  for(i in approved){
    udv <- svd(LS[[effs[i]]])
    scores[[effs[i]]]    <- (udv$u * rep(udv$d, each=N))[,1:sum(assign==i), drop=FALSE]
    dimnames(scores[[effs[i]]]) <- list(rownames(LS[[effs[i]]]), paste("Comp", 1:sum(assign==i), sep=" "))
    loadings[[effs[i]]]  <- udv$v[,1:sum(assign==i), drop=FALSE]
    dimnames(loadings[[effs[i]]]) <- list(colnames(LS[[effs[i]]]), paste("Comp", 1:sum(assign==i), sep=" "))
    projected[[effs[i]]] <- residuals %*% loadings[[effs[i]]]
    dimnames(projected[[effs[i]]]) <- list(rownames(LS[[effs[i]]]), paste("Comp", 1:sum(assign==i), sep=" "))
    singulars[[effs[i]]] <- udv$d[1:sum(assign==i)]
    names(singulars[[effs[i]]]) <- paste("Comp", 1:sum(assign==i), sep=" ")
    if(pca.in!=0){ # Transform back if PCA on Y has been performed
      loadings[[effs[i]]] <- Yudv$v[,1:pca.in,drop=FALSE] %*% loadings[[effs[i]]]
      dimnames(loadings[[effs[i]]]) <- list(colnames(LS[[effs[i]]]), paste("Comp", 1:sum(assign==i), sep=" "))
    }
  }
  
  # Reset options
  options(opt)
  obj <- list(scores=scores, loadings=loadings, projected=projected, singulars=singulars, 
              LS=LS, effects=effects, Y=Y, X=M, residuals=residuals,
              ssq=unlist(ssq), ssqY=ssqY, explvar=unlist(ssq)/ssqY,
              call=match.call(), fit.type=fit.type)
  if(pca.in!=0){
    obj$Ypca <- list(svd=Yudv, ncomp=pca.in)
  }
  class(obj) <- c('asca', 'list')
  return(obj)
}

#' @rdname asca
#' @export
print.asca <- function(x, ...){
  cat("Anova Simultaneous Component Analysis fitted using", x$fit.type)
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  invisible(x)
}
  
#' @rdname asca
#' @export
summary.asca <- function(object, ...){
  warning('Not implemented yet!')
}

#' @rdname asca
#' @export
loadings.asca <- function(object, factor = 1, ...){
  loads <- object$loadings[[factor]]
  class(loads) <- "loadings"
  return(loads)
}
#' @rdname asca
#' @export
loadingplot.asca <- function(object, factor = 1, comps = 1:2, ...){
  plot(loadings(object=object, factor=factor), comps=comps, ...)
}

#' @rdname asca
#' @export
scores.asca <- function(object, factor = 1, ...){
  scors <- object$scores[[factor]]
  class(scors) <- "scores"
  return(scors)
}

#' @rdname asca
#' @export
projections <- function (object, ...) {
  UseMethod("projections", object)
}
#' @rdname asca
#' @export
projections.asca <- function(object, factor = 1, ...){
  projs <- object$projected[[factor]]
  class(projs) <- "projs"
  return(projs)
}

#' @rdname asca
#' @export
scoreplot.asca <- function(object, factor = 1, comps = 1:2, pch.scores = 19, pch.projections = 1, 
                           gr.col = 1:nlevels(object$effects[[factor]]), ellipsoids,
                           xlim,ylim, xlab,ylab, legendpos, ...){
  # Number of levels in current factor
  nlev  <- nlevels(object$effects[[factor]])
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
      points(scors[as.numeric(object$effects[[factor]]) == i, comps], pch=pch.scores, col=gr.col[i])
      points(projs[as.numeric(object$effects[[factor]]) == i, comps], pch=pch.projections, col=gr.col[i])
    }
    if(!missing(legendpos))
      legend(legendpos, legend = levels(object$effects[[factor]]), col=gr.col, pch=pch.scores)

    if(!missing(ellipsoids)){
      if(ellipsoids == "data"){
        dataEllipse(projs[,comps], groups = object$effects[[factor]], levels=c(0.4,0.68,0.95), add=TRUE, plot.points=FALSE, col=gr.col, lwd=1, group.labels="", center.pch=FALSE, lty=3)
        # for(i in 1:nlev){
        #   sigma <- cov(projs[as.numeric(object$effects[[factor]]) == i, comps])
        #   S     <- chol(sigma)
        #   nobji <- sum(as.numeric(object$effects[[factor]]) == i)
        #   
        # }
      }
    }
  } else { # Line plot
    plot(scors[,comps], as.numeric(object$effects[[factor]]), xlim=xlim, 
         ylim=ylim, xlab=xlab, ylab=ylab, axes = FALSE)
    axis(1)
    axis(2, at=1:nlev, labels = levels(object$effects[[factor]]))
    box()
    for(i in 1:nlev){
      points(scors[as.numeric(object$effects[[factor]]) == i, comps], rep(i,sum(as.numeric(object$effects[[factor]]) == i)), pch=pch.scores, col=gr.col[i])
      points(projs[as.numeric(object$effects[[factor]]) == i, comps], rep(i,sum(as.numeric(object$effects[[factor]]) == i)), pch=pch.projections, col=gr.col[i])
    }
  }
  #TODO: 
  # - Data ellipsoids, car::dataEllipse
  # - Model ellipsoids
  
}
