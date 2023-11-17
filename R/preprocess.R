#' @name preprocess
#' @aliases block.preprocess preprocess
#' @title Preprocessing of block data
#' 
#' @description
#' This is an interface to simplify preprocessing of one, a subset or all
#' blocks in a multiblock object, e.g., a \code{data.frame} (see \code{block.data.frame})
#' or \code{list}. Several standard preprocessing methods are supplied in addition to 
#' letting the user supply it's own function.
#'
#' @param X \code{data.frame} or \code{list} of data.
#' @param block \code{vector} of block(s) to preprocess (\code{integer}s or \code{character}s).
#' @param fun \code{character} or \code{function} selecting which preprocessing to apply (see Details).
#' @param ... additional arguments to underlying functions.
#' 
#' @details
#' The \code{fun} parameter controls the type of preprocessing to be performed:
#' * autoscale: centre and scale each feature/variable.
#' * center: centre each feature/variable.
#' * scale: scale each feature/variable.
#' * SNV: Standard Normal Variate correction, i.e., centre and scale each sample across features/variables.
#' * EMSC: Extended Multiplicative Signal Correction defaulting to basic EMSC (2nd order polynomials). Further parameters are sent to \code{EMSC::EMSC}.
#' * Fro: Frobenius norm scaling of whole block.
#' * FroSq: Squared Frobenius norm scaling of whole block (sum of squared values).
#' * SingVal: Singular value scaling of whole block (first singular value).
#' * User defined: If a function is supplied, this will be applied to chosen blocks. 
#' Preprocessing can be done for all blocks or a subset. It can also be done in a series of operations to combine preprocessing techniques.
#' 
#' @return The input multiblock object is preprocessed and returned.
#'
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results and plotting are found in \code{\link{multiblock_results}} and \code{\link{multiblock_plots}}, respectively.
#' 
#' @examples
#' data(potato)
#' # Autoscale Chemical block
#' potato <- block.preprocess(potato, block = "Chemical", "autoscale")
#' # Apply SNV to NIR blocks
#' potato <- block.preprocess(potato, block = 3:4, "SNV")
#' # Centre Sensory block
#' potato <- block.preprocess(potato, block = "Sensory", "center")
#' # Scale all blocks to unit Frobenius norm
#' potato <- block.preprocess(potato, fun = "Fro")
#' 
#' # Effect of SNV
#' NIR <- (potato$NIRraw + rnorm(26)) * rnorm(26,1,0.2)
#' NIRc <- block.preprocess(list(NIR), fun = "SNV")[[1]]
#' old.par <- par(mfrow = c(2,1), mar = c(4,4,1,1))
#' matplot(t(NIR), type="l", main = "uncorrected", ylab = "")
#' matplot(t(NIRc), type="l", main = "corrected", ylab = "")
#' par(old.par)
#' 
#' @export
block.preprocess <- function(X, block = 1:length(X), 
                             fun = c("autoscale", "center", "scale", "SNV", 
                                      "EMSC", "Fro", "FroSq", "SingVal"), ...){
  if(inherits(fun, "function")){
    user.fun <- fun
    fun <- "user"
  }
    
  fun <- fun[1]
  if(!(fun %in% c("autoscale", "center", "scale", "SNV", 
                   "EMSC", "Fro", "FroSq", "SingVal", "user")))
     stop("Unknown preprocessing chosen")
  for(b in block){
    if(fun == "autoscale")
      X[[b]] <- scale(X[[b]])
    if(fun == "center")
      X[[b]] <- scale(X[[b]], scale = FALSE)
    if(fun == "scale")
      X[[b]] <- scale(X[[b]], center = FALSE)
    if(fun == "SNV")
      X[[b]] <- t(scale(t(X[[b]])))
    if(fun == "EMSC"){
        if("EMSC" %in% rownames(installed.packages())){
          X[[b]] <- EMSC::EMSC(X[[b]], ...)$corrected
        } else {
        stop("To preprocess with 'EMSC', please install the 'EMSC' package, e.g., using\ninstall.packages('EMSC')\nand rerun this function\n")
      }
    }
    if(fun == "Fro")
      X[[b]] <- X[[b]]/norm(X[[b]], 'F')
    if(fun == "Fro2")
      X[[b]] <- X[[b]]/norm(X[[b]], 'F')^2
    if(fun == "SingVal")
      X[[b]] <- X[[b]]/svd(X[[b]])$d[1]
    if(fun == "user")
      X[[b]] <- user.fun(X[[b]])
  }
  return(X)
}