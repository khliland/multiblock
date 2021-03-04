#' @name asca
#' @aliases asca
#' @title Analysis of Variance Simultaneous Component Analysis - ASCA
#'
#' @param formula Model formula accepting a single response (block) and predictor names separated by + signs.
#' @param data The data set to analyse.
#' @param subset Expression for subsetting the data before modelling.
#' @param weights Optional object weights.
#' @param subset Subset of objects
#' @param na.action How to handle NAs (no action implemented).
#' @param family Error distributions and link function for Generalized Linear Models.
#' @param pca.in Compress response before ASCA (number of components).
#'
#' @return An \code{asca} object containing loadings, scores, explained variances, etc.
#' 
#' @description ASCA is a method which decomposes a multivariate response according to one or more design
#' variables. ANOVA is used to split variation into contributions from factors, and PCA is performed
#' on the corresponding least squares estimates, i.e., $Y = X1 B1 + X2 B2 + ... + E = T1 P1' + T2 P2' + ... + E$.
#' 
#' @references 
#' * Smilde, A., Jansen, J., Hoefsloot, H., Lamers,R., Van Der Greef, J., and Timmerman, M.(2005). ANOVA-Simultaneous Component Analysis (ASCA): A new tool for analyzing designed metabolomics data. Bioinformatics, 21(13), 3043–3048.
#' * Liland, K.H., Smilde, A., Marini, F., and Næs,T. (2018). Confidence ellipsoids for ASCA models based on multivariate regression theory. Journal of Chemometrics, 32(e2990), 1–13.
#' * Martin, M. and Govaerts, B. (2020). LiMM-PCA: Combining ASCA+ and linear mixed models to analyse high-dimensional designed data. Journal of Chemometrics, 34(6), e3232.
#'
#' @importFrom lme4 lmer
#' @importFrom car ellipse dataEllipse
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' @examples
#' # Simulate data set
#' dataset   <- data.frame(y = I(matrix(rnorm(24*10),ncol=10)), 
#'                         x = factor(c(rep(2,8),rep(1,8),rep(0,8))), 
#'                         z = factor(rep(c(1,0),12)), 
#'                         w = rnorm(24))
#' colnames(dataset$y) <- paste('Var', 1:10, sep=" ")
#' rownames(dataset) <- paste('Obj', 1:24, sep=" ")
#' 
#' # Basic ASCA model with two factors
#' mod <- asca(y~x+z, data=dataset)
#' print(mod)
#' 
#' # Result plotting for first factor
#' loadingplot(mod, scatter=TRUE)
#' scoreplot(mod)
#' 
#' # ASCA model with compressed response using 5 principal components
#' mod.pca <- asca(y~x+z, data=dataset, pca.in=5)
#' 
#' # Mixed Model ASCA
#' mod <- asca(y~x+(1|z), data=dataset)
#' 
#' @export
asca <- function(formula, data, subset, weights, na.action, family, pca.in = FALSE){
  ## Force contrast to sum
  opt <- options(contrasts = c(unordered="contr.sum", ordered="contr.poly"))
  
  ## Get the data matrices
  Y <- data[[formula[[2]]]]
  N <- nrow(Y)
  p <- ncol(Y)
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
    maxDir <- min(sum(assign==i), p)
    udv <- svd(LS[[effs[i]]])
    scores[[effs[i]]]    <- (udv$u * rep(udv$d, each=N))[,1:maxDir, drop=FALSE]
    dimnames(scores[[effs[i]]]) <- list(rownames(LS[[effs[i]]]), paste("Comp", 1:maxDir, sep=" "))
    loadings[[effs[i]]]  <- udv$v[,1:maxDir, drop=FALSE]
    dimnames(loadings[[effs[i]]]) <- list(colnames(LS[[effs[i]]]), paste("Comp", 1:maxDir, sep=" "))
    projected[[effs[i]]] <- residuals %*% loadings[[effs[i]]]
    dimnames(projected[[effs[i]]]) <- list(rownames(LS[[effs[i]]]), paste("Comp", 1:maxDir, sep=" "))
    singulars[[effs[i]]] <- udv$d[1:maxDir]
    names(singulars[[effs[i]]]) <- paste("Comp", 1:maxDir, sep=" ")
    if(pca.in!=0){ # Transform back if PCA on Y has been performed
      loadings[[effs[i]]] <- Yudv$v[,1:pca.in,drop=FALSE] %*% loadings[[effs[i]]]
      dimnames(loadings[[effs[i]]]) <- list(colnames(LS[[effs[i]]]), paste("Comp", 1:maxDir, sep=" "))
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
  # # Experimental features
  # # Generalised ASCA, here with a mock Gaussian distribution
  # mod.glm <- asca(y~x+z, data=dataset, family="gaussian")
  # 
  # # Generalised Mixed Model ASCA
  # mod <- asca(y~x+(1|z), data=dataset, family="gaussian")
}
