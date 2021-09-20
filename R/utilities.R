# Create all unique component combinations of n_blocks from min_level to max_level in length
#' Unique combinations of blocks
#' 
#' @description Compute a list of all possible block combinations where
#' the number of blocks in each combination is limited by parameters
#' \code{min_level} and \code{max_level}.
#'
#' @param n_block \code{integer} number of input blocks.
#' @param max_level \code{integer} maximum number of blocks per combination.
#' @param min_level \code{integer} minimum number of blocks per combination.
#' 
#' @details This function is used for minimal redundancy implementations of
#' \code{\link{rosa}} and \code{\link{sopls}} and for lookups into computed
#' components.
#'
#' @return A list of unique block combinations.
#'
#' @examples
#' unique_combos(3, 2)
#' 
#' @export
unique_combos <- function(n_block, max_level, min_level=2){
  all_comb_lim <- function(comb, vec, nb, n, pos){
    if(pos > n){
      comb <- c(comb, sort(vec))
    } else {
      for(i in 1:nb){
        vec[pos] <- i
        comb <- all_comb_lim(comb, vec, nb, n, pos+1)
      }
    }
    return(comb)
  }
  comb <- unique(matrix(all_comb_lim(numeric(0), numeric(0), n_block, max_level, 1), ncol=max_level, byrow = TRUE))
  comb <- unique(apply(comb,1,unique))
  comb <- comb[order(unlist(lapply(comb,length)), decreasing = TRUE)]
  if(min_level>1){
    comb <- comb[unlist(lapply(comb,length))>=min_level]
  }
  return(comb)
}

pcaopt <- function(X,T,P,ncomp){
  ResVar <- numeric(ncomp+1)
  ResVar[1] <- mean(X^2)
  for(i in 1:ncomp){
    E <- X - tcrossprod(T[,1:i,drop=FALSE], P[,1:i,drop=FALSE])
    ResVar[i+1] <- mean(E^2)
  }
  crit <- (1:ncomp)*0.02*ResVar[1] + ResVar[-1]
  which.min(crit)
}

normCols <- function(X){
  X <- as.matrix(X)
  norms <- sqrt(colSums(X*X))
  X / rep(norms, each=nrow(X))
}

normVec <- function(x){
  x / sqrt(sum(crossprod(x)))
}

colnamesList <- function(X, nam){
  for(i in 1:length(X))
    colnames(X[[i]]) <- nam
  return(X)
}

rownamesList <- function(X, nam){
  for(i in 1:length(X))
    rownames(X[[i]]) <- nam
  return(X)
}


#' Colour palette generation from matrix of RGB values
#'
#' @param n Integer number of colorus to produce.
#' @param colmatrix A numeric \code{matrix} of three columns (R,G,B) to generate colour palette from.
#'
#' @return A vector of n colours in hexadecimal RGB format.
#'
#' @examples 
#' mcolors(5)
#' @export
mcolors <- function(n, colmatrix = matrix(c(0,0,1, 1,1,1, 1,0,0), 3,3, byrow = TRUE)){
  cols <- character(n)
  colR <- approx(seq(0,1,length.out=nrow(colmatrix)), colmatrix[,1], (0:(n-1))/(n-1))$y
  colG <- approx(seq(0,1,length.out=nrow(colmatrix)), colmatrix[,2], (0:(n-1))/(n-1))$y
  colB <- approx(seq(0,1,length.out=nrow(colmatrix)), colmatrix[,3], (0:(n-1))/(n-1))$y
  for(i in 1:n){
    cols[i] <- rgb(colR[i],colG[i],colB[i])
  }
  cols
}


whichMins <- function(short, long){
  out <- integer(length(short))
  for(i in 1:length(short)){
    out[i] <- which.min((short[i]-long)^2)
  }
  out[short < min(long)] <- NA; out[short > max(long)] <- NA
  out
}

#' @title Block-wise indexable data.frame
#' 
#' @description This is a convenience function for making \code{data.frame}s that are easily
#' indexed on a block-wise basis.
#'
#' @param X Either a single \code{data.frame} to index or a \code{list} of matrices/data.frames
#' @param block_inds Named \code{list} of indexes if \code{X} is a single \code{data.frame}, otherwise \code{NULL}.
#' @param to.matrix \code{logical} indicating if input list elements should be converted to matrices.
#'
#' @return A \code{data.frame} which can be indexed block-wise.
#' @examples
#' # Random data
#' M <- matrix(rnorm(200), nrow = 10)
#' # .. with dimnames
#' dimnames(M) <- list(LETTERS[1:10], as.character(1:20))
#' 
#' # A named list for indexing
#' inds <- list(B1 = 1:10, B2 = 11:20)
#' 
#' X <- block.data.frame(M, inds)
#' str(X)
#' 
#' @export
block.data.frame <- function(X, block_inds = NULL, to.matrix = TRUE){
  if(!is.null(block_inds)){
    # Use indices to convert matrix/data.frame into list of matrices/data.frames
    if(is.null(names(block_inds)))
      warning("When 'block_inds' is supplied, it should be a named list with indices/names of variables associated with blocks.")
    Z <- lapply(block_inds, function(i)X[,i,drop=FALSE])
    X <- lapply(Z, function(z){rownames(z) <- rownames(X);z})
  }
  # Enclose blocks in "as.is"
  if(to.matrix)
    X <- lapply(X, function(x)I(as.matrix(x)))
  else
    X <- lapply(X, function(x)I(x))
  # Return as data.frame
  X <- do.call(data.frame, X)
  X <- lapply(X, function(x){rownames(x) <- rownames(X);x})
  return(X)
}


#' Dummy-coding of a single vector
#' 
#' @description Flexible dummy-coding allowing for all R's built-in types of contrasts
#' and optional dropping of a factor level to reduce rank defficiency probability.
#'
#' @param Y \code{vector} to dummy code.
#' @param contrast Contrast type, default = "contr.sum".
#' @param drop \code{logical} indicating if one level should be dropped (default = TRUE).
#'
#' @return \code{matrix} made by dummy-coding the input vector.
#'
#' @examples
#' vec <- c("a","a","b","b","c","c")
#' dummycode(vec)
#' @export
dummycode <- function(Y, contrast = "contr.sum", drop = TRUE){
  nlev <- nlevels(Y)
  lev  <- levels(Y)
  if(drop){
    X    <- model.matrix(~x,data.frame(x=factor(Y)), contrasts.arg = list(x=contrast))
    X    <- X[, -1, drop=FALSE]
  } else {
    X    <- model.matrix(~x-1,data.frame(x=factor(Y)), contrasts.arg = list(x=contrast))
  }
  attributes(X) <- list(dim = attributes(X)$dim, dimnames = attributes(X)$dimnames)
  X
}


#' Explained predictor variance
#' 
#' @description Extraction and/or computation of explained variances for various
#' object classes in the \code{multiblock} package.
#'
#' @param object An object fitted using a method from the multiblock package
#'
#' @return A \code{vector} of component-wise explained variances for predictors.
#' @examples
#' data(potato)
#' so <- sopls(Sensory ~ Chemical + Compression, data=potato, ncomp=c(10,10), 
#'             max_comps=10)
#' explvar(so)
#' @export
explvar <- function(object){
  switch(class(object)[1],
         rosa = 100 * object$Xvar / object$Xtotvar,
         sopls = apply(crossprod(x<-scale(do.call(cbind,object$data$X), scale=FALSE),object$decomp$T)^2,2,sum)/sum(x^2)*100,
         multiblock = object$explvar,
         mvr = 100 * object$Xvar / object$Xtotvar,
         princomp =,
         prcomp = 100 * object$sdev^2 / sum(object$sdev^2),
         scores = attr(object, "explvar"),
         loadings = attr(object, "explvar"),
         scores.multiblock = attr(object, "explvar")
  )}

#' Vector of component names
#'
#' @description Convenience function for creating a vector
#' of component names based on the dimensions the input object
#' (\code{matrix} or object having a \code{score} function).
#'
#' @param object An object fitted using the multiblock package.
#' @param comps \code{integer} vector of components.
#' @param explvar \code{logical} indicating if explained variances should be included.
#' @param ... Unused
#'
#' @return A \code{character} vector of component names.
#' 
#' @details This is a copy of \code{compnames} from the \code{pls} package to work with
#' \code{multiblock} objects.
#' @export
compnames <- function(object, comps, explvar = FALSE, ...) {
  M <- if (is.matrix(object)) object else scores(object)
  labs <- colnames(M)
  if (missing(comps))
    comps <- seq(along = labs)
  else
    labs <- labs[comps]
  if (isTRUE(explvar) && !is.null(evar <- explvar(M)[comps]))
    labs <- paste(labs, " (", format(evar, digits = 2, trim = TRUE),
                  " %)", sep = "")
  return(labs)
}