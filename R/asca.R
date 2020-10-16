
#' Title
#'
#' @param formula 
#' @param data 
#' @param subset 
#' @param weights 
#' @param na.action 
#' @param family 
#' @param pca.in 
#'
#' @return
#'
#' @examples
#' dataset1   <- data.frame(y=I(matrix(rnorm(24*10),ncol=10)), x=factor(c(rep(2,8),rep(1,8),rep(0,8))), z=factor(rep(c(1,0),12)), w=rnorm(24))
#' mod <- asca(y~x+z, data=dataset1)
#' mod <- asca(y~x+z, data=dataset1, pca.in=0)
#' mod <- asca(y~x+z+w, data=dataset1, family="gaussian")
#' mod <- asca(y~x+(1|z), data=dataset1)
#' mod <- asca(y~x+(1|z)+w, data=dataset1, family="gaussian")
#' @export
asca <- function(formula, data, subset, weights, na.action, family, pca.in = FALSE){
  ## Force contrast to sum
  opt <- options(contrasts = c(unordered="contr.sum", ordered="contr.poly"))
  
  ## Get the data matrices
  Y <- data[[formula[[2]]]]
  N <- nrow(Y)
  Y <- Y - rep(colMeans(Y), each=N) # Centre Y
  if(pca.in != 0){
    if(pca.in == 1)
      stop('pca.in = 1 is not supported (single response)')
    Yudv <- svd(Y)
    Y <- Yudv$u[,1:pca.in,drop=FALSE] * rep(Yudv$d[1:pca.in], each=N)
  }
  residuals <- Y
  
  mf <- match.call(expand.dots = FALSE)
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
    }
  }
  M      <- model.matrix(ano)
  effs   <- attr(terms(ano), "term.labels")
  assign <- attr(M, "assign")
  
  # Exclude numeric effects and their interactions
  nums   <- names(which(unlist(lapply(model.frame(ano), class)) == "numeric"))
  if(length(nums)>0){
    exclude  <- match(nums, rownames(attr(terms(ano), "factors")))
    approved <- which(colSums(attr(terms(ano), "factors")[exclude,,drop=FALSE])==0)
  } else {
    approved <- 1:max(assign)
  }
  LS <- list()
  for(i in approved){
    LS[[effs[i]]] <- M[, assign==i, drop=FALSE] %*% coefs[assign==i,]
    residuals <- residuals - LS[[effs[i]]]
  }
  
  # SCAs
  scores <- loadings <- projected <- singulars <- list()
  for(i in approved){
    udv <- svd(LS[[effs[i]]])
    scores[[effs[i]]]    <- (udv$u * rep(udv$d, each=N))[,1:sum(assign==i), drop=FALSE]
    loadings[[effs[i]]]  <- udv$v[,1:sum(assign==i), drop=FALSE]
    projected[[effs[i]]] <- residuals %*% loadings[[effs[i]]]
    singulars[[effs[i]]] <- udv$d[1:sum(assign==i)]
    if(pca.in!=0) # Transform back if PCA on Y has been performed
      loadings[[effs[i]]] <- Yudv$v[,1:pca.in,drop=FALSE] %*% loadings[[effs[i]]]
  }
  
  # Reset options
  options(opt)
  obj <- list(scores=scores, loadings=loadings, projected=projected, 
              singulars=singulars, LS=LS, Y=Y, X=M, residuals=residuals)
  if(pca.in!=0){
    obj$Ypca <- list(svd=Yudv, ncomp=pca.in)
  }
  class(obj) <- c('ASCA', 'list')
  return(obj)
}
