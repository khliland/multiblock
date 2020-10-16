# Create all unique component combinations of n_blocks from min_level to max_level in length
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


#' Color palette generation from matrix of RGB values
#'
#' @param n Integer number of colors to produce.
#' @param colmatrix A numeric \code{matrix} of three columns (R,G,B) to generate color palette from.
#'
#' @return A vector of n colors.
#' @export
#'
#' @examples
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
