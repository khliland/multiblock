#' L-PLS data simulation for exo-type analysis
#'
#' @param I \code{numeric} number of rows of X1 and X2
#' @param N \code{numeric} number of columns in X1 and X3
#' @param J \code{numeric} number of columns in X2
#' @param K \code{numeric} number of rows in X3
#' @param ncomp \code{numeric} number of latent components
#' @author Solve Sæbø (adapted by Kristian Hovde Liland)
#' 
#' @description Three data blocks are simulated to express covariance in an exo-L-PLS direction (see \code{\link{lpls}}.
#' Dimensionality and number of underlying components can be controlled.
#'
#' @return A \code{list} of three matrices with dimensions matching in an L-shape.
#' 
#' @examples
#' lp <- lplsData(I = 30, N = 20, J = 5, K = 6, ncomp = 2)
#' names(lp)
#' 
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' @export
lplsData <- function(I = 30, N = 20, J = 5, K = 6, ncomp = 2){
  # Data-simulation for LPLS testing (exo-type)
  
  # Simulations
  X1rowm <- rnorm(I, 0, 1)
  X1colm <- rnorm(N, 2, 2)
  X2colm <- rnorm(J, 5, 1)
  X3colm <- rnorm(N, 0, 1)
  X3rowm <- rnorm(K, 2, 1)
  
  X1 <- matrix(rnorm(I*N,0,1),I,N)
  ULV <- svd(X1)
  T21 <- ULV$u[,1:ncomp]%*%diag(sqrt(ULV$d[1:ncomp]))
  T22 <- ULV$v[,1:ncomp]%*%diag(sqrt(ULV$d[1:ncomp]))
  P1 <- matrix(rnorm(ncomp*J,4,1),J,ncomp)
  P3 <- matrix(rnorm(ncomp*K,2,2),K,ncomp)
  X1 <- T21%*%t(T22)
  X2 <- T21%*%t(P1)
  X3 <- T22%*%t(P3)
  
  #   X3<-t(X3)
  
  dimnames(X2) <- list(paste("I",1:I,sep=""), paste("class",1:J))
  dimnames(X1) <- list(paste("I",1:I,sep=""), as.character(1:N))
  dimnames(X3) <- list(as.character(1:N), paste("clust",1:K))
  
  return(list(X1 = X1, X2 = X2, X3 = t(X3)))
}
