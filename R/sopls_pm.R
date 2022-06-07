########################################################
# Total, direct, indirect, additional, overall effects #
########################################################
# ... between first and last block given intermediates

#' @name SO_TDI
#' @aliases sopls_pm print.SO_TDI
#' @title Total, direct, indirect and additional effects in SO-PLS-PM.
#'
#' @param X A \code{list} of input blocks (of type \code{matrix}).
#' @param Y A \code{matrix} of response(s).
#' @param ncomp An \code{integer} vector giving the number of components per block or a single integer for common number of components.
#' @param max_comps Maximum total number of components.
#' @param sel.comp A \code{character} or \code{integer} vector indicating the type ("opt" - minimum error / "chi" - chi-squared reduced) or exact number of components in selections.
#' @param computeAdditional A \code{logical} indicating if additional components should be computed.
#' @param sequential A \code{logical} indicating if sequential component optimization should be applied.
#' @param B An \code{integer} giving the number of bootstrap replicates for variation estimation.
#' @param k An \code{integer} indicating number of cross-validation segments (default = 10).
#' @param type A \code{character} for selecting type of cross-validation segments (default = "consecutive").
#' @param simultaneous \code{logical} indicating if simultaneous orthogonalisation on intermediate blocks should be performed (default = TRUE).
#' @param x An object of type \code{SO_TDI}.
#' @param showComp A \code{logical} indicating if components should be shown in print (default = TRUE).
#' @param heading A \code{character} giving the heading of the print.
#' @param digits An \code{integer} for selecting number of digits in print.
#' @param ... Not implemented
#' 
#' @description SO-PLS-PM is the use of SO-PLS for path-modelling. This particular function
#' is used to compute effects (explained variances) in sub-paths of the directed acyclic graph.
#' 
#' @details \code{sopls_pm} computes 'total', 'direct', 'indirect' and 'additional' effects for the 'first' versus the 
#' 'last' input block by cross-validated explained variances. 'total' is the explained variance when doing
#' regression of 'first' -> 'last'. 'indirect' is the the same, but controlled for the intermediate blocks.
#' 'direct' is the left-over part of the 'total' explained variance when subtracting the 'indirect'. Finally,
#' 'additional' is the added explained variance of 'last' for each block following 'first'.
#' 
#' \code{sopls_pm_multiple} is a wrapper for \code{sopls_pm} that repeats the calculation for all pairs of blocks
#' from 'first' to 'last'. Where \code{sopls_pm} has a separate response, Y, signifying the 'last' block, 
#' \code{sopls_pm_multiple} has multiple 'last' blocks, depending on sub-path, thus collects the response(s)
#' from the list of blocks X.
#' 
#' When sel.comp = "opt", the number of components for all models are optimized using cross-validation
#' within the ncomp and max_comps supplied. If sel.comp is "chi", an optimization is also performed,
#' but parsimonious solutions are sought through a chi-square chriterion. When setting sel.comp to a
#' numeric vector, exact selection of number of components is performed.
#' 
#' When setting B to a number, e.g. 200, the procedures above are repeated B times using bootstrapping
#' to estimate standard deviations of the cross-validated explained variances.
#'
#' @return An object of type \code{SO_TDI} containing total, direct and indirect effects, plus
#' possibly additional effects and standard deviations (estimated by bootstrapping).
#' 
#' @references 
#' * Menichelli, E., Almøy, T., Tomic, O., Olsen, N. V., & Næs, T. (2014). SO-PLS as an exploratory tool for path modelling. Food quality and preference, 36, 122-134.
#' * Næs, T., Romano, R., Tomic, O., Måge, I., Smilde, A., & Liland, K. H. (2020). Sequential and orthogonalized PLS (SO-PLS) regression for path analysis: Order of blocks and relations between effects. Journal of Chemometrics, e3243.
#'
#' @examples
#' # Single path for the potato data:
#' data(potato)
#' pot.pm <- sopls_pm(potato[1:3], potato[['Sensory']], c(5,5,5), computeAdditional=TRUE)
#' pot.pm
#' 
#' # Corresponding SO-PLS model:
#' # so <- sopls(Sensory ~ ., data=potato[c(1,2,3,9)], ncomp=c(5,5,5), validation="CV", segments=10)
#' # maageSeq(pot.so, compSeq = c(3,2,4))
#' 
#' # All path in the forward direction for the wine data:
#' data(wine)
#' pot.pm.multiple <- sopls_pm_multiple(wine, ncomp = c(4,2,9,8))
#' pot.pm.multiple
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' @export
sopls_pm <- function(X, Y, ncomp, max_comps = min(sum(ncomp), 20), sel.comp = "opt", computeAdditional = FALSE, sequential = FALSE, B = NULL, k = 10, type = "consecutive", simultaneous = TRUE){#validation = "LOO", ...){
  Y <- as.matrix(Y)
  n      <- dim(Y)[1]
  nresp  <- dim(Y)[2]
  nblock <- length(X)
  
  ###############
  # Check inputs
  if(is.numeric(sel.comp) && sum(sel.comp)>max_comps){
    stop(paste0("'sel.comp' is outside bounds set by 'max_comps'. Suggested max_comp = ", sum(sel.comp)))
  }
  ## Check components
  if(length(ncomp)==1 && length(X)>1){
    ncomp <- rep(ncomp, length(X))
  }
  comps <- pmin(ncomp, max_comps)
  if(is.numeric(sel.comp) && any(sel.comp>comps)){
    stop("'sel.comp' is larger than 'comps'")
  }
  if(nblock==1){
    computeAdditional <- FALSE
  }
  
  ##############
  # CV segments
  segments <- cvsegments(N=n, k=k, type=type)
  
  ###############
  # TOTAL effect
  mod_total    <- plsr(Y~X[[1]], ncomp = min(comps[1],max_comps), validation="CV", segments=segments)
  r <- matrix(RMSEP(mod_total)$val[1,,], nresp)
  # r[,1] <- MSEP0(Y,segments)
  denom <- apply(Y,2,var)*(n-1)/n
  R2_total_ind <- 1-r^2/denom # r[,1]^2             # PRESS0 variant
  # R2_total_ind <- drop(R2(mod_total)$val) # SST variant
  R2_total <- 1 - colSums(r^2)/sum(denom)
  # R2_total     <- apply(R2_total_ind,2,mean)
  if(is.character(sel.comp)){
    if(sel.comp == "opt"){
      comp_total <- which.max(R2_total)-1
      total <- max(R2_total)
      total_ind <- R2_total_ind[,comp_total+1]
    } else { # "chi"
      r <- sqrt(colMeans(r^2))
      comp_total <- chi2cv(r,n,0.05)-1
      total <- R2_total[comp_total+1]
      total_ind <- R2_total_ind[,comp_total+1]
    }
  } else {
    comp_total <- sel.comp[1]
    total <- R2_total[sel.comp[1]+1]
    total_ind <- R2_total_ind[,sel.comp[1]+1]
  }
  
  # Bootstrap total effect
  if(!is.null(B)){
    Boots <- matrix(0L, n, B)
    for(i in 1:B){
      Boots[,i] <- sample.int(n)
    }
    total_boot <- bootPLSR(X[[1]],Y, max(comp_total,1), Boots, segments=segments)
  }
  
  ################
  # DIRECT effect
  first <- X[[1]] - rep(colMeans(X[[1]]), each=n)
  last <- Y - rep(colMeans(Y), each=n)
  if(nblock>1){
    S_first <- S_last <- list()
    for(i in 2:nblock){
      # First block
      if(comps[i]>0){
        mod_first <- plsr(first~X[[i]], ncomp = min(comps[i],max_comps), validation="CV", segments=segments)
        r      <- matrix(RMSEP(mod_first)$val[1,,], ncol(first))
        # r[,1]  <- MSEP0(first,segments)
        denom  <- apply(first,2,var)*(n-1)/n
        r2_ind <- 1-r^2/denom # r[,1]^2
        r2     <- 1 - colSums(r^2)/sum(denom)
        # r2     <- apply(r2_ind,2,mean)
        
        if(is.character(sel.comp) && sel.comp == "chi"){
          r <- sqrt(colMeans(r^2))
          ind <- chi2cv(r,n,0.05)
        } else { # "opt" or numeric
          ind <- which.max(r2)
        }
      } else {
        ind <- 1
      }
      if(ind > 1){
        S <- mod_first$scores[,1:(ind-1),drop=FALSE]
        S_first[[i-1]] <- S
        if(!simultaneous)
          first <- first - S %*% tcrossprod(solve(crossprod(S)), S) %*% first
      }
      
      # Last block
      if(comps[i]>0){
        mod_last <- plsr(last~X[[i]], ncomp = min(comps[i],max_comps), validation="CV", segments=segments)
        r      <- matrix(RMSEP(mod_last)$val[1,,], nresp)
        # r[,1]  <- MSEP0(last,segments)
        denom  <- apply(last,2,var)*(n-1)/n
        r2_ind <- 1-r^2/denom # r[,1]^2
        r2     <- 1 - colSums(r^2)/sum(denom)
        # r2     <- apply(r2_ind,2,mean)
        if(is.character(sel.comp) && sel.comp == "chi"){
          r <- sqrt(colMeans(r^2))
          ind <- chi2cv(r,n,0.05)
        } else { # "opt" or numeric
          ind <- which.max(r2)
        }
      } else {
        ind <- 1
      }
      if(ind > 1){
        S <- mod_last$scores[,1:(ind-1),drop=FALSE]
        S_last[[i-1]] <- S
        if(!simultaneous)
          last <- last - S %*% tcrossprod(solve(crossprod(S)), S) %*% last
      }
    }
    S_first <- do.call(cbind,S_first)
    S_last <- do.call(cbind,S_last)
    if(simultaneous){
      if(!is.null(S_first))
        first <- first - S_first %*% tcrossprod(solve(crossprod(S_first)), S_first) %*% first
      if(!is.null(S_last))
        last <- last - S_last %*% tcrossprod(solve(crossprod(S_last)), S_last) %*% last
    }
  }
  # Last ~ First
  mod_first_last <- plsr(last~first, ncomp = min(comps[1],max_comps), validation="CV", segments=segments)
  r      <- matrix(RMSEP(mod_first_last)$val[1,,], nresp)
  # r[,1]  <- MSEP0(last,segments)
  denom <- apply(last,2,var)*(n-1)/n
  r2_ind <- 1-r^2/denom # r[,1]^2
  r2     <- 1 - colSums(r^2)/sum(denom)
  # r2     <- apply(r2_ind,2,mean)
  rescale<-  mean((last-rep(colMeans(last), each=n))^2) / mean((Y-rep(colMeans(Y), each=n))^2)
  if(is.character(sel.comp) && sel.comp == "chi"){
    r <- sqrt(colMeans(r^2))
    comp_direct <- chi2cv(r,n,0.05)-1
    direct <- r2[comp_direct+1] * rescale
  } else { # "opt" or numeric
    direct <- max(r2) * rescale
    comp_direct <- which.max(r2)-1
  }
  direct_ind <- r2_ind[, comp_direct+1] * rescale
  direct <- max(min(total, direct),0)
  direct_ind <- pmax(pmin(total_ind, direct_ind),0)
  
  # Bootstrap direct effect
  if(!is.null(B)){
    direct_boot <- bootPLSR(first,last, max(comp_direct,1), Boots, segments=segments)*rescale
  }
  
  ##################
  # INDIRECT effect
  indirect <- total - direct
  indirect_ind <- total_ind - direct_ind
  
  if(computeAdditional){
    #################
    # OVERALL effect
    if(sequential){
      mod_overall <- sopls_cv_seq(X=X, Y=Y, comps=comps, max_comps = max_comps, sel.comp = sel.comp, segments=segments)
      R2_overall <- mod_overall$expl_var
      R2_overall_ind <- mod_overall$expl_var_ind
      comp_overall <- match(paste0(mod_overall$chosen,collapse=","), names(mod_overall$RMSECV))
      overall <- R2_overall[comp_overall]
      overall_ind <- R2_overall_ind[,comp_overall]
    } else {
      mod_overall <- sopls_cv(X = X, Y = Y, comps = comps, max_comps = max_comps, segments=segments)
      R2_overall <- mod_overall$expl_var
      R2_overall_ind <- mod_overall$expl_var_ind
      if(is.character(sel.comp)){
        if(sel.comp == "opt"){
          comp_overall <- which.max(R2_overall)
          overall <- max(R2_overall)
          overall_ind <- R2_overall_ind[,comp_overall]
        } else { # "chi"
          comp_overall <- chi2cv(mod_overall$RMSECV,n,0.05)
          overall <- R2_overall[comp_overall]
          overall_ind <- R2_overall_ind[,comp_overall]
        }
      } else {
        comp_overall <- match(paste0(sel.comp,collapse=","), names(R2_overall))
        overall <- R2_overall[comp_overall]
        overall_ind <- R2_overall_ind[,comp_overall]
      }
    }
    the_comps <- as.numeric(strsplit(names(R2_overall)[comp_overall],",")[[1]])
    if(!is.null(B)){
      overall_boot <- bootSOPLS(X,Y, the_comps, max_comps, Boots, segments=segments)
    }
    
    ####################
    # ADDITIONAL effect
    comp_seq <- the_comps # as.numeric(strsplit(names(overall),",")[[1]])
    additional <- numeric(nblock)
    mod_add <- mod_overall #sopls_cv(X, Y, comps=comp_seq, max_comps=max_comps, validation=validation, ...)
    curr <- rep(0, nblock)
    for(i in 1:nblock){
      curr[i] <- comp_seq[i]
      ind <- match(paste(curr, collapse=","), names(mod_add$RMSECV))
      additional[i] <- mod_add$expl_var[ind]
    }
    additional <- diff(additional)
    if(!is.null(B)){
      additional_boot <- overall_boot # bootSOPLS(X,Y, comp_seq, max_comps, Boots, segments=segments)
      additional_boot <- apply(additional_boot,2,sd)
    }
  }
  
  ##################
  # Returned object
  out <- list(effects = list(direct=unname(direct), indirect=unname(indirect), total=unname(total)),
              models = list(total=mod_total), comps = list(total=comp_total, direct=comp_direct))
  if(computeAdditional){
    out$effects$additional <- unname(additional)
    out$effects$overall    <- unname(overall)
    out$models$overall     <- mod_overall 
    out$comp$additional    <- comp_seq
  }
  if(!is.null(B)){
    out$sd <- list(direct=sd(direct_boot), indirect=sd(total_boot-direct_boot), total=sd(total_boot))
    if(computeAdditional){
      out$sd$additional <- additional_boot
      out$sd$overall    <- sd(overall_boot)
    }
  }
  
  class(out) <- list("SO_TDI", "list")
  out
}


#' @rdname SO_TDI
#' @export
print.SO_TDI <- function(x, showComp=TRUE, heading="SO-PLS path effects", digits=2, ...){
  # cat(paste0(heading,"\n"))
  if(!is.null(x$sd)){
    comp_text <- paste(",", unlist(x$comp[1:2]))
    dfr <- data.frame(direct = paste(round(x$effects$direct*100, digits=digits), " (", round(x$sd$direct*100, digits=digits),ifelse(showComp,comp_text[2],""),")", sep="" ),
                     indirect = paste(round(x$effects$indirect*100, digits=digits), " (", round(x$sd$indirect*100, digits=digits),")", sep="" ),
                     total = paste(round(x$effects$total*100, digits=digits), " (", round(x$sd$total*100, digits=digits),ifelse(showComp,comp_text[1],""),")", sep="" ))
    if(length(x$effects$additional)>0){
      comp_text <- paste(",", x$comp[[3]])
      for(i in 1:length(x$effects$additional)){
        dfr <- cbind(dfr, additional = paste(round(x$effects$additional[i]*100, digits=digits), " (", round(x$sd$additional[i]*100, digits=digits),ifelse(showComp,comp_text[i],""),")", sep="" ))
      }
      dfr <- cbind(dfr, overall = paste(round(x$effects$overall*100, digits=digits), " (", round(x$sd$total*100, digits=digits),")", sep="" ))
      cn <- colnames(dfr)
      if(length(x$effects$additional)>1){
        cn[4:(length(cn)-1)] <- paste0("additional", 1:length(x$effects$additional))
      } else {
        cn[4:(length(cn)-1)] <- "additional"
      }
      colnames(dfr) <- cn
    }
    rownames(dfr) <- ""
  } else {
    comp_text <- paste(" (", unlist(x$comp[1:2]), ")", sep="")
    dfr <- data.frame(direct = paste(round(x$effects$direct*100, digits=digits),ifelse(showComp,comp_text[2],""), sep="" ),
                      indirect = paste(round(x$effects$indirect*100, digits=digits), sep="" ),
                      total = paste(round(x$effects$total*100, digits=digits),ifelse(showComp,comp_text[1],""), sep="" ))
    if(length(x$effects$additional)>0){
      comp_text <- paste(" (", x$comp[[3]], ")", sep="")
      for(i in 1:length(x$effects$additional)){
        dfr <- cbind(dfr, additional = paste(round(x$effects$additional[i]*100, digits=digits),ifelse(showComp,comp_text[i],""), sep="" ))
      }
      dfr <- cbind(dfr, overall = paste(round(x$effects$overall*100, digits=digits), sep="" ))
      cn <- colnames(dfr)
      if(length(x$effects$additional)>1){
        cn[4:(length(cn)-1)] <- paste0("additional", 1:length(x$effects$additional))
      } else {
        cn[4:(length(cn)-1)] <- "additional"
      }
      colnames(dfr) <- cn
    }
    rownames(dfr) <- ""
  }
  print(dfr)
}


#############################################
# Multiple calls for each sequential version
#' @rdname SO_TDI
#' @export
sopls_pm_multiple <- function(X, ncomp, max_comps = min(sum(ncomp), 20), sel.comp = "opt", computeAdditional = FALSE, sequential = FALSE, B = NULL, k = 10, type = "consecutive"){
  nblock <- length(X)
  
  ## Check components
  if(length(ncomp)==1 && length(X)>1){
    ncomp <- rep(ncomp, length(X))
  }
  comps <- pmin(ncomp, max_comps)

  if(is.null(X_names <- names(X)))
    X_names <- paste0('Block ', 1:nblock)
  models <- list()
  model_names <- list()
  b <- 1
  for(first in 1:(nblock-1)){
    for(last in (first+1):nblock){
      if(is.numeric(sel.comp)){
        sc <- sel.comp[first:(last-1)]
      } else {
        sc <- sel.comp
      }
      models[[b]] <- sopls_pm(X[first:(last-1)], X[[last]], comps[first:(last-1)], max_comps = max_comps, sel.comp = sc, computeAdditional = computeAdditional, sequential = sequential, B = B, k = k, type = type)
      model_names[[b]] <- paste0(X_names[first],"->",X_names[last])
      b <- b+1
    }
  }
  names(models) <- model_names
  class(models) <- list("SO_TDI_multiple", "list")
  models
}

#' @rdname SO_TDI
#' @export
print.SO_TDI_multiple <- function(x, heading="SO-PLS path effects", digits=2, ...){
#  cat(paste0(heading,"\n"))
  print(x[])
}
#   if(!is.null(x[[1]]$direct_sd)){
#     m <- as.data.frame(matrix(unlist(x), ncol=10, byrow=TRUE))
#     colnames(m) <- names(x[[1]])
#     rownames(m) <- names(x)
#     dfr <- data.frame(direct = paste(round(m$direct[1]*100, digits=digits), " (", round(m$direct_sd[1]*100, digits=digits),")", sep="" ),
#                       indirect = paste(round(m$indirect[1]*100, digits=digits), " (", round(m$indirect_sd[1]*100, digits=digits),")", sep="" ),
#                       total = paste(round(m$total[1]*100, digits=digits), " (", round(m$total_sd[1]*100, digits=digits),")", sep="" ),
#                       additional = paste(round(m$additional[1]*100, digits=digits), " (", round(m$additional_sd[1]*100, digits=digits),")", sep="" ),
#                       overall = paste(round(m$overall[1]*100, digits=digits), " (", round(m$overall_sd[1]*100, digits=digits),")", sep="" ))
#     for(i in 2:dim(m)[1]){
#       dfr <- rbind(dfr, data.frame(direct = paste(round(m$direct[i]*100, digits=digits), " (", round(m$direct_sd[i]*100, digits=digits),")", sep="" ),
#                                    indirect = paste(round(m$indirect[i]*100, digits=digits), " (", round(m$indirect_sd[i]*100, digits=digits),")", sep="" ),
#                                    total = paste(round(m$total[i]*100, digits=digits), " (", round(m$total_sd[i]*100, digits=digits),")", sep="" ),
#                                    additional = paste(round(m$additional[i]*100, digits=digits), " (", round(m$additional_sd[i]*100, digits=digits),")", sep="" ),
#                                    overall = paste(round(m$overall[i]*100, digits=digits), " (", round(m$overall_sd[i]*100, digits=digits),")", sep="" )))
#     }
#     rownames(dfr) <- rownames(m)
#   }
#   print(dfr)
# }


####################################
# Model-based bootstrapping of PLSR
bootPLSR <- function(X,Y, ncomp, B, segments){
  n   <- dim(Y)[1]
  nB  <- dim(B)[2]
  R2B <- numeric(nB)
  nresp <- dim(Y)[2]
  
  mod_base <- plsr(Y ~ X, ncomp=ncomp, validation="CV", segments=segments)
  fits <- drop(predict(mod_base, ncomp=ncomp))
  res  <- Y-fits
  for(b in 1:nB){
    Yb <- fits + res[B[,b],,drop=FALSE]
    mod_boot <- plsr(Yb ~ X, ncomp=ncomp, validation="CV", segments=segments)
    r <- matrix(RMSEP(mod_boot)$val[1,,], nresp)
    # r[,1] <- MSEP0(Yb,segments)
    denom <- apply(Yb,2,var)*(n-1)/n
    r2 <- 1-r^2/denom # r[,1]^2
    R2B[b] <- max(apply(r2,2,mean))
  }
  R2B
}

######################################
# Model-based bootstrapping of SO-PLS
bootSOPLS <- function(X,Y, comps, max_comps, B, segments){
  n   <- dim(Y)[1]
  nB  <- dim(B)[2]
  nresp  <- dim(Y)[2]
  nblock <- length(X)
  R2B    <- matrix(0,nB,nblock)

  mod_base <- sopls_single_prediction(X, X, Y, comps, max_comps)
  fits <- mod_base$Ypred
  fits <- fits[,,dim(fits)[3]]
  res  <- Y-fits
  for(b in 1:nB){
    Yb <- fits + res[B[,b],,drop=FALSE]
    mod_boot <- sopls_cv(X, Yb, comps=comps, max_comps=max_comps, segments=segments)
    curr <- rep(0, nblock)
    for(i in 1:nblock){
      curr[i] <- comps[i]
      ind <- match(paste(curr, collapse=","), names(mod_boot$RMSECV))
      R2B[b,i] <- mod_boot$expl_var[ind]
    }
  }
  R2B
}

#########
# MSEP0 
MSEP0 <- function(Y, segments){
  Y  <- as.matrix(Y)
  np <- dim(Y)
  Yhat <- matrix(0.0, np[1], np[2])
  for(i in 1:length(segments)){
    Yhat[segments[[i]], ] <- matrix(colMeans(Yhat[-segments[[i]], ]), nrow=length(segments[[i]]), ncol=np[2], byrow=TRUE)
  }
  return(colMeans((Y-Yhat)^2))
}


###################
# Function testing
# DITAO <- sopls_pm(X, Y, c(3,2,2,1), max_comps = 8, sel.comp = "opt", sequential = FALSE, validation = "CV", k=10)
# print(DITAO)
# DITAOchi <- sopls_pm(X, Y, c(3,2,2,1), max_comps = 8, sel.comp = "chi", sequential = FALSE, validation = "CV", k=10)
# print(" ");print(DITAOchi)
# DITAOseq <- sopls_pm(X, Y, c(3,2,2,1), max_comps = 8, sel.comp = "opt", sequential = TRUE, validation = "CV", k=10)
# print(" ");print(DITAOseq)
# DITAOseqchi <- sopls_pm(X, Y, c(3,2,2,1), max_comps = 8, sel.comp = "chi", sequential = TRUE, validation = "CV", k=10)
# print(" ");print(DITAOseqchi)
# 
# 
# DITAOB <- sopls_pm(X, Y, c(3,2,2,1), max_comps = 8, sel.comp = "opt", sequential = FALSE, B=10, validation = "CV", k=10)
# print(DITAOB)
# 
# ########################
# 
# 
# 
# ADE <- sopls_pm(process_data[c(1,4)], process_data[[5]], c(4,3), sel.comp = c(4,3), sequential = TRUE, B=10)
# BDE <- sopls_pm(process_data[c(2,4)], process_data[[5]], c(12,4), sel.comp = c(12,4), sequential = TRUE, B=10)
# CDE <- sopls_pm(process_data[c(3,4)], process_data[[5]], c(9,2), sel.comp = c(9,2), sequential = TRUE, B=10)
# 
# ADE1 <- sopls_pm(process_data[c(1,4)], process_data[[5]][,1], c(3,3), sel.comp = "opt", sequential = TRUE, B=10)
# BDE1 <- sopls_pm(process_data[c(2,4)], process_data[[5]][,1], c(5,3), sel.comp = "opt", sequential = TRUE, B=10)
# CDE1 <- sopls_pm(process_data[c(3,4)], process_data[[5]][,1], c(3,1), sel.comp = "opt", sequential = TRUE, B=10)
# 
# ADE2 <- sopls_pm(process_data[c(1,4)], process_data[[5]][,2], c(3,3), sel.comp = "opt", sequential = TRUE, B=10)
# BDE2 <- sopls_pm(process_data[c(2,4)], process_data[[5]][,2], c(5,3), sel.comp = "opt", sequential = TRUE, B=10)
# CDE2 <- sopls_pm(process_data[c(3,4)], process_data[[5]][,2], c(3,1), sel.comp = "opt", sequential = TRUE, B=10)
# 
# ADE3 <- sopls_pm(process_data[c(1,4)], process_data[[5]][,3], c(3,3), sel.comp = "opt", sequential = TRUE, B=10)
# BDE3 <- sopls_pm(process_data[c(2,4)], process_data[[5]][,3], c(5,3), sel.comp = "opt", sequential = TRUE, B=10)
# CDE3 <- sopls_pm(process_data[c(3,4)], process_data[[5]][,3], c(3,1), sel.comp = "opt", sequential = TRUE, B=10)
# 
# print(ADE, "A+D->E")
# print(BDE, "B+D->E")
# print(CDE, "C+D->E")
# print(ADE1, "A+D->E1")
# print(BDE1, "B+D->E1")
# print(CDE1, "C+D->E1")
# print(ADE2, "A+D->E2")
# print(BDE2, "B+D->E2")
# print(CDE2, "C+D->E2")
# print(ADE3, "A+D->E3")
# print(BDE3, "B+D->E3")
# print(CDE3, "C+D->E3")


