#' @keywords internal
#' @noRd
notin <- function(a,b){
  a[!(a %in% b)]
}

#' @keywords internal
#' @noRd
lmR2 <- function(formula, data, ...){
  data <- as.data.frame(data)
  N <- nrow(data[[1]])
  Nseq <- seq(1,N)
  pred <- as.matrix(predict(lm(formula, data, ...)))
  Y <- model.response(model.frame(formula,data))
  SSE <- colSums((Y-pred)^2)
  SST <- colSums((scale(Y,scale=FALSE)^2))
  return(list(R2 = 1-sum(SSE)/sum(SST), R2ind = 1-(SSE/SST)))
}

#' @keywords internal
#' @noRd
lmR2cv <- function(formula, data, segments = 5, type = "consecutive", ...){
  data <- as.data.frame(data)
  N <- nrow(data[[1]])
  Nseq <- seq(1,N)
  pred <- 0
  cv <- pls::cvsegments(N = N, k = segments, type = type, ...)
  for(i in 1:segments){
    train <- match(Nseq, cv[[i]], nomatch=0) == 0
    modi <- lm(formula, data[train, ], ...)
    if(length(pred) == 1){
      pred1 <- predict(modi, newdata=data[cv[[i]],])
      pred <- matrix(0, nrow = N, ncol = ifelse(is.vector(pred1), 1, ncol(pred1)))
      pred[cv[[i]],] <- pred1
    } else {
      pred[cv[[i]],] <- predict(modi, newdata=data[cv[[i]],])
    }
  }
  Y <- model.response(model.frame(formula,data))
  SSE <- colSums((Y-pred)^2)
  SST <- colSums((scale(Y,scale=FALSE)^2))
  return(list(R2 = 1-sum(SSE)/sum(SST), R2ind = 1-(SSE/SST)))
}

#' @keywords internal
#' @noRd
#' @importFrom pls plsr R2
#' @importFrom pracma Rank
TotUnCoAd_regression <- function(relations, start, end, blocks, ncomp, validation, segments, SO=TRUE, fits=FALSE, ...){
  # Convert input dictionary to lists
  b_names <- names(blocks)
  block_list <- unname(blocks)
  
  # Convert start/end to index if needed
  if (is.character(start)) {
    start <- which(b_names == start)
  }
  if (is.character(end)) {
    end <- which(b_names == end)
  }
  
  # First and last blocks
  start_block <- as.matrix(block_list[[start]])
  end_block <- as.matrix(block_list[[end]])
  
  # Find direct paths to end (R equivalent of find_incomming)
  find_incoming <- function(relations, target) {
    relations[,1][apply(relations, 1, function(row) target == row[2])]
  }
  
  to_end <- find_incoming(relations, end)
  to_end_no_start <- to_end[to_end != start]

  # ncomp semantics:
  # - positive values: maximum number of components per predictor block
  # - negative values: force number of components per predictor block
  # ncomp is provided as a per-block vector (non-predictor blocks may be 0)
  has_ncomp <- !missing(ncomp) && !is.null(ncomp)
  force_mode <- FALSE
  start_req <- others_req <- all_req <- 0
  if(has_ncomp){
    force_mode <- any(ncomp < 0)
    start_raw <- if(start <= length(ncomp)) ncomp[start] else 0
    others_raw <- if(length(to_end_no_start) > 0) sum(ncomp[to_end_no_start], na.rm = TRUE) else 0
    all_raw <- start_raw + others_raw
    start_req <- abs(start_raw)
    others_req <- abs(others_raw)
    all_req <- abs(all_raw)
  }
  
  if (length(to_end_no_start) > 0) {
    others <- as.matrix(do.call(cbind, block_list[to_end_no_start]))
    #alls <- do.call(cbind, block_list[to_end])
    alls <- cbind(start_block, others)
    noOthers <- FALSE
  } else {
    alls <- as.matrix(start_block)
    #alls <- do.call(cbind, block_list[to_end])
    noOthers <- TRUE
  }
  
  # All
  df <- data.frame(all = I(alls),
                   end = I(end_block))
  if(SO){
    ncompa_max <- min(nrow(alls)-1, ncol(alls))
    ncompa <- ncompa_max
    if(has_ncomp && all_req > 0){
      if(force_mode){
        if(all_req > ncompa_max){
          stop("Forced ncomp for 'all' exceeds feasible maximum for start=", b_names[start], " and end=", b_names[end])
        }
        ncompa <- all_req
      } else {
        ncompa <- min(all_req, ncompa_max)
      }
    }
    pls <- pls::plsr(end ~ all, data = df, ncomp = ncompa, 
                     validation = validation, segments = segments, ...)
    ncompa <- dim(pls$fitted.values)[3]
    SSEcv <- cbind("(Intercept)"=pls$validation$PRESS0, pls$validation$PRESS)
    SST <- colSums(scale(df$end,scale=FALSE)^2)
    ev <- 1 - colSums(SSEcv)/sum(SST)
    if(has_ncomp && force_mode && all_req > 0){
      ncomp_all <- all_req + 1
    } else {
      ncomp_all <- which.max(ev)
    }
    if(fits){
      r2 <- c("(Intercept)"=pls::R2(pls, estimate="train")$val[1], 
              1-apply((rep(df$end,dim(pls$fitted.values)[3])-pls$fitted.values)^2,3,sum)/
                sum((df$end-rep(colMeans(df$end),each=nrow(df$end)))^2))
      SSE <- SSEcv*0
      for(i in 1:ncompa)
        SSE[,i] <- colSums((df$end - pls$fitted.values[,,i])^2)
      All_EF <- max(0, r2[names(ncomp_all)])
      if(r2[names(ncomp_all)]<=0){
        ncomp_all <- 0
        All_EF_ind <- rep(0, ncol(df$end))
      } else {
        All_EF_ind <- pmax(0, (1-SSE/SST)[, ncomp_all])
      }
    } else {
      All_EF <- max(0,max(ev))
      if(max(ev)<=0){
        ncomp_all <- 0
        All_EF_ind <- rep(0, ncol(df$end))
      } else {
        All_EF_ind <- pmax(0, (1-SSEcv/SST)[, ncomp_all])
      }
    }    
  } else {
    if(fits){
      reg <- lmR2(end ~ all, data = df)
    } else {
      reg <- lmR2cv(end ~ all, data = df, segments = segments, ...)
    }
    All_EF <- reg$R2
    All_EF_ind <- reg$R2ind
  }
  
  # Only unique
  if(noOthers){
    # Relation without alternative paths
    df <- data.frame(start = I(start_block),
                     end = I(end_block))
    if(SO){
      ncomp_unique_max <- min(nrow(start_block)-1, ncol(start_block))
      ncomp_unique_fit <- ncomp_unique_max
      if(has_ncomp && start_req > 0){
        if(force_mode){
          if(start_req > ncomp_unique_max){
            stop("Forced ncomp for 'start' exceeds feasible maximum for start=", b_names[start], " and end=", b_names[end])
          }
          ncomp_unique_fit <- start_req
        } else {
          ncomp_unique_fit <- min(start_req, ncomp_unique_max)
        }
      }
      pls <- pls::plsr(end ~ start, data = df, ncomp = ncomp_unique_fit, 
                       validation = validation, segments = segments, ...)
      ncomp_unique_fit <- dim(pls$fitted.values)[3]
      SSEcv <- cbind("(Intercept)"=pls$validation$PRESS0, pls$validation$PRESS)
      SST <- colSums(scale(df$end,scale=FALSE)^2)
      ev <- 1 - colSums(SSEcv)/sum(SST)
      if(has_ncomp && force_mode && start_req > 0){
        ncomp_unique <- start_req + 1
      } else {
        ncomp_unique <- which.max(ev)
      }
      if(fits){
        r2 <- c("(Intercept)"=pls::R2(pls, estimate="train")$val[1], 
                1-apply((rep(df$end,dim(pls$fitted.values)[3])-pls$fitted.values)^2,3,sum)/
                  sum((df$end-rep(colMeans(df$end),each=nrow(df$end)))^2))
        SSE <- SSEcv*0
        for(i in 1:ncomp_unique_fit)
          SSE[,i] <- colSums((df$end - pls$fitted.values[,,i])^2)
        Tot_EF <- max(0, r2[names(ncomp_unique)])
        r2_ind <- pls::R2(pls, estimate = "train")$val
        if(r2[names(ncomp_unique)]<=0){
          ncomp_unique <- 0
          Tot_EF_ind <- rep(0, ncol(df$end))
        } else {
          Tot_EF_ind <- pmax(0, r2_ind[,,ncomp_unique])
        }
      } else {
        Tot_EF <- max(0,max(ev))
        if(max(ev)<=0){
          ncomp_unique <- 0
          Tot_EF_ind <- rep(0, ncol(df$end))
        } else {
          Tot_EF_ind <- pmax(0, (1-SSEcv/SST)[, ncomp_unique])
        }
      }
    } else {
      if(fits){
        reg <- lmR2(end ~ start, data = df)
      } else {
        reg <- lmR2cv(end ~ start, data = df, segments = segments, ...)
      }
      Tot_EF <- max(reg$R2, 0)
      Tot_EF_ind <- pmax(reg$R2ind,0)
    }
    miss_U <- FALSE
    if(!any(relations[,1]==start & relations[,2]==end)){
      miss_U <- TRUE
    }
    ret <- list(Tot_EF = Tot_EF, Un_EF = Tot_EF, Co_EF = NaN, Ad_EF = NaN, other_EF = NaN, All_EF = All_EF,
                Tot_EF_ind = Tot_EF_ind, Un_EF_ind = Tot_EF_ind, Co_EF_ind = NaN, Ad_EF_ind = NaN, 
                other_EF_ind = NaN, All_EF_ind = All_EF_ind,
                other_names = list(), miss_U=miss_U, vars = colnames(df$end))
    if(SO){
      ret$ncomps <- c(ncomp_unique = ncomp_unique, ncomp_all = ncomp_all)-1
    }
    class(ret) <- c("TotUnCoAd_regression", class(ret))
    return(ret)
    
  } else {
    # Data for models
    df <- data.frame(start = I(start_block),
                     others = I(others),
                     all = I(cbind(start_block, others)),
                     end = I(end_block))
    if(SO){
      seg_size <- ceiling(nrow(start_block)/segments)
      ncomp_start_max <- min(nrow(start_block)-1-seg_size, pracma::Rank(start_block))
      ncomp_others_max <- min(nrow(others)-1-seg_size, pracma::Rank(others))
      if(has_ncomp && (start_req > 0 || others_req > 0)){
        req_start <- if(start_req > 0) start_req else ncomp_start_max
        req_others <- if(others_req > 0) others_req else ncomp_others_max
        if(force_mode){
          if(req_start > ncomp_start_max){
            stop("Forced ncomp for 'start' exceeds feasible maximum for start=", b_names[start], " and end=", b_names[end])
          }
          if(req_others > ncomp_others_max){
            stop("Forced ncomp for 'others' exceeds feasible maximum for start=", b_names[start], " and end=", b_names[end])
          }
          ncomp_pair <- c(req_start, req_others)
        } else {
          ncomp_pair <- c(min(req_start, ncomp_start_max), min(req_others, ncomp_others_max))
        }
      } else {
        ncomp_pair <- c(ncomp_start_max, ncomp_others_max)
      }
      # SO-PLS for Tot_EF and Ad_EF
      so_T_Ad <- sopls(end ~ start + others, data=df, ncomp = ncomp_pair, max_comps = sum(ncomp_pair),
                       validation = validation, segments=segments, sequential=TRUE, ...)
      # Explained variance
      if(force_mode){
        chosen <- ncomp_pair
      } else {
        chosen <- so_T_Ad$validation$chosen
      }
      ev <- so_T_Ad$validation$expl_var
      ev1 <- ev[grep(",0", names(ev))]
      ncomp_total <- chosen[1]
      names(ncomp_total) <- paste(chosen[1],",0", sep="")
      ev2 <- ev[which(unlist(lapply(strsplit(names(ev), ","), 
                                    function(i)i[1]))==strsplit(names(ncomp_total),",")[[1]][1])]
      ncomp_additional <- chosen[2]
      names(ncomp_additional) <- paste(chosen[1],",",chosen[2], sep="")
      if(fits){
        r2 <- pls::R2(so_T_Ad, estimate = "train", individual=FALSE)
        r2_ind <- pls::R2(so_T_Ad, estimate = "train", individual=TRUE)
      } else {
        r2 <- pls::R2(so_T_Ad, estimate = "CV", individual=FALSE)
        r2_ind <- pls::R2(so_T_Ad, estimate = "CV", individual=TRUE)
      }
      if(names(ncomp_total) == "0,0"){
        ev1 <- so_T_Ad$validation$expl_var[1]
        if(fits){
          ev1 <- ev1*0
        }
      } else {
        ev1 <- r2[names(ncomp_total)]
      }
      if(names(ncomp_additional) == "0,0"){
        ev2 <- so_T_Ad$validation$expl_var[1]
        if(fits){
          ev2 <- ev2*0
        }
      } else {
        ev2 <- r2[names(ncomp_additional)]
      }
      if(is.matrix(r2_ind)){
        if(names(ncomp_total) == "0,0"){
          ev1_ind <- so_T_Ad$validation$expl_var_ind[,1]
          if(fits){
            ev1_ind <- ev1_ind*0
          }
        } else {
          ev1_ind <- r2_ind[,names(ncomp_total)]
        }
        if(names(ncomp_additional) == "0,0"){
          ev2_ind <- so_T_Ad$validation$expl_var_ind[,1]
          if(fits){
            ev2_ind <- ev2_ind*0
          }
        } else {
          ev2_ind <- r2_ind[,names(ncomp_additional)]
        }
      } else {
        if(names(ncomp_total) == "0,0"){
          ev1_ind <- so_T_Ad$validation$expl_var_ind[1]
          if(fits){
            ev1_ind <- ev1_ind*0
          }
        } else {
          ev1_ind <- r2_ind[names(ncomp_total)]
        }
        if(names(ncomp_additional) == "0,0"){
          ev2_ind <- so_T_Ad$validation$expl_var_ind[1]
          if(fits){
            ev2_ind <- ev2_ind*0
          }
        } else {
          ev2_ind <- r2_ind[names(ncomp_additional)]
        }
      }
      Tot_EF <- pmax(0,ev1)
      Tot_EF_ind <- pmax(0,ev1_ind)
      Ad_EF <- pmax(0,ev2-ev1)
      Ad_EF_ind <- pmax(0,ev2_ind-ev1_ind)
      
      # SO-PLS for Other_EF and Un_EF
      ncomp_rev <- ncomp_pair[2:1]
      so_Others_U <- sopls(end ~ others + start, data=df, ncomp = ncomp_rev, max_comps = sum(ncomp_rev),
                           validation = validation, segments=segments, sequential=TRUE, ...)
      # Explained variance
      if(force_mode){
        chosen <- ncomp_rev
      } else {
        chosen <- so_Others_U$validation$chosen
      }
      ev <- so_Others_U$validation$expl_var
      ev1 <- ev[grep(",0", names(ev))]
      ncomp_others <- chosen[1]
      names(ncomp_others) <- paste(chosen[1],",0", sep="")
      ev2 <- ev[which(unlist(lapply(strsplit(names(ev), ","), 
                                    function(i)i[1]))==strsplit(names(ncomp_others),",")[[1]][1])]
      ncomp_unique <- chosen[2]
      names(ncomp_unique) <- paste(chosen[1],",",chosen[2], sep="")
      if(fits){
        r2 <- pls::R2(so_Others_U, estimate = "train", individual=FALSE)
        r2_ind <- pls::R2(so_Others_U, estimate = "train", individual=TRUE)
      } else {
        r2 <- pls::R2(so_Others_U, estimate = "CV", individual=FALSE)
        r2_ind <- pls::R2(so_Others_U, estimate = "CV", individual=TRUE)
      }
      if(names(ncomp_others) == "0,0"){
        ev1 <- so_Others_U$validation$expl_var[1]
        if(fits){
          ev1 <- ev1*0
        }
      } else {
        ev1 <- r2[names(ncomp_others)]
      }
      if(names(ncomp_unique) == "0,0"){
        ev2 <- so_Others_U$validation$expl_var[1]
        if(fits){
          ev2 <- ev2*0
        }
      } else {
        ev2 <- r2[names(ncomp_unique)]
      }
      if(is.matrix(r2_ind)){
        if(names(ncomp_others) == "0,0"){
          ev1_ind <- so_Others_U$validation$expl_var_ind[,1]
          if(fits){
            ev1_ind <- ev1_ind*0
          }
        } else {
          ev1_ind <- r2_ind[,names(ncomp_others)]
        }
        if(names(ncomp_unique) == "0,0"){
          ev2_ind <- so_Others_U$validation$expl_var_ind[,1]
          if(fits){
            ev2_ind <- ev2_ind*0
          }
        } else {
          ev2_ind <- r2_ind[,names(ncomp_unique)]
        }
      } else {
        if(names(ncomp_others) == "0,0"){
          ev1_ind <- so_Others_U$validation$expl_var_ind[1]
          if(fits){
            ev1_ind <- ev1_ind*0
          }
        } else {
          ev1_ind <- r2_ind[names(ncomp_others)]
        }
        if(names(ncomp_unique) == "0,0"){
          ev2_ind <- so_Others_U$validation$expl_var_ind[1]
          if(fits){
            ev2_ind <- ev2_ind*0
          }
        } else {
          ev2_ind <- r2_ind[names(ncomp_unique)]
        }
      }
      other_EF <- pmax(0,ev1)
      other_EF_ind <- pmax(0,ev1_ind)
      Un_EF <- pmax(0,ev2-ev1)
      Un_EF_ind <- pmax(0,ev2_ind-ev1_ind)
    } else { # OLS
      # Total effect
      if(fits){
        reg <- lmR2(end ~ start, data = df)
      } else {
        reg <- lmR2cv(end ~ start, data = df, segments = segments, ...)
      }
      Tot_EF <- reg$R2
      Tot_EF_ind <- reg$R2ind
      # Effect of others
      if(fits){
        reg <- lmR2(end ~ others, data = df)
      } else {
        reg <- lmR2cv(end ~ others, data = df, segments = segments, ...)
      }
      other_EF <- reg$R2
      other_EF_ind <- reg$R2ind
      Un_EF <- max(0, All_EF - other_EF)
      Un_EF_ind <- pmax(0, All_EF_ind - other_EF_ind)
      Ad_EF <- max(0, All_EF - Tot_EF)
      Ad_EF_ind <- pmax(0, All_EF_ind - Tot_EF_ind)
    }
    if(SO){
      Co_EF <- Tot_EF - Un_EF
      Co_EF_ind <- Tot_EF_ind - Un_EF_ind
    } else {
      Co_EF <- All_EF - Un_EF - Ad_EF
      Co_EF_ind <- All_EF_ind - Un_EF_ind - Ad_EF_ind
    }
    
    miss_U <- FALSE
    if(!any(relations[,1]==start & relations[,2]==end)){
      miss_U <- TRUE
    }
    ret <- list(Tot_EF = Tot_EF, Un_EF = Un_EF, Co_EF = Co_EF, Ad_EF = Ad_EF, other_EF = other_EF, All_EF = All_EF,
                Tot_EF_ind = Tot_EF_ind, Un_EF_ind = Un_EF_ind, Co_EF_ind = Co_EF_ind, Ad_EF_ind = Ad_EF_ind, 
                other_EF_ind = other_EF_ind, All_EF_ind = All_EF_ind, 
                other_names = b_names[to_end_no_start], miss_U = miss_U, vars = colnames(df$end))
    if(SO){
      ret$ncomps <- c("ncomp_all" = ncomp_all-1, "ncomp_total" = ncomp_total, 
                      "ncomp_additional" = ncomp_additional, "ncomp_others" = ncomp_others, 
                      "ncomp_unique" = ncomp_unique)
    }
    class(ret) <- c("TotUnCoAd_regression", class(ret))
    return(ret)
  }
}

#' Analyse Multivariate Path Effects
#'
#' Decomposes effects between blocks in a multiblock system into total effects, unique effects,
#' common contributions, and additional effects. Supports both SO-PLS (Sequential Orthogonalised PLS)
#' and ordinary least squares (OLS) approaches with cross-validation or fitted values.
#'
#' @param relations A matrix defining the path structure. Rows specify relationships
#'   between blocks, with 2 columns: `c(from_block_index, to_block_index)`.
#' @param blocks A multiblock data frame or list of named matrices containing the block data.
#' @param validation Validation method passed to [pls::plsr()]. Common options:
#'   `"CV"` for cross-validation, "LOO" for leave-one-out validation, or `"none"` for fitted values.
#' @param segments Number of segments for cross-validation. Default is 5.
#' @param SO Logical. If `TRUE` (default), use SO-PLS; if `FALSE`, use OLS.
#' @param fits Logical. If `TRUE`, use fitted values instead of cross-validation.
#'   Default is `FALSE`.
#' @param boot Number of bootstrap samples for computing standard errors.
#'   Default is 0 (no bootstrapping).
#' @param ncomp Optional component settings per predictor block.
#'   Supply either:
#'   - a vector with length equal to the number of predictor blocks (`unique(relations[,1])`),
#'     in ascending predictor-block order, or
#'   - a full-length vector with one entry per block (non-predictors may be `0`).
#'   Positive integers set upper limits; negative integers force exact component counts.
#'   All predictor entries must be non-zero and have the same sign.
#' @param transitive_closure Logical. If `TRUE`, automatically add transitive closure to the relations. Default is `FALSE`.
#' @param ... Additional arguments passed to underlying fitting functions.
#'
#' @return An object of class `path_effects`, which is a matrix with the following
#'   components for each path:
#'   - `Tot`: Total effect
#'   - `Un`: Unique effect  
#'   - `Co`: Common contribution
#'   - `Ad`: Additional effect
#'
#'   Attributes include:
#'   - `scaled`: Scaled contribution matrix
#'   - `individual`: Individual-level contributions
#'   - `boot`: Bootstrap replicates (if `boot > 0`)
#'
#' @seealso [print.path_effects()], [plot.path_effects()]
#'
#' @examples
#' # Analysis of the mobile dataset
#' data(mobile)
#' 
#' # Define path structure (A->B, A->E, A->G, B->C, B->D, B->E, C->D, D->E, 
#' #                        D->E, E->F, E->G, F->G)
#' paths <- matrix(c(1,2, 1,5, 1,7, 2,3, 2,4, 2,5, 3,4, 3,5, # Add 0,2, for A->C
#'                   4,5, 5,6, 5,7, 6,7),
#'                 ncol=2, byrow=TRUE)
#' 
#' # Compute path effects with cross-validation using SO-PLS
#' pem <- path_effects(paths, mobile, validation="CV", segments=5, 
#'                     segment.type="consecutive")
#' 
#' # Print results
#' print(pem)
#' 
#' # Plot all results
#' plot(pem)
#' 
#' # Print and plot single path
#' print(pem, "A","G")
#' plot(pem, from = "A", to = "G")
#'
#' # Print and plot results per variable
#' print(pem, individual = TRUE)
#' plot(pem, individual = TRUE)
#' 
#' # Analysis of the NIR-Raman-PUFA data (emulsions) 
#' data(emulsions)
#' 
#' # Standardise response
#' emulsions$PUFA <- scale(emulsions$PUFA)
#' 
#' # Define path structure (NIR->Raman, NIR->PUFA, Raman->PUFA)
#' paths_NRP <- matrix(c(1,2,1,3,2,3), ncol = 2, byrow = TRUE)
#' 
#' \dontrun{ # Too time consuming
#' # Compute path effects with cross-validation using SO-PLS
#' pem_NRP <- path_effects(paths_NRP, emulsions, validation="CV", 
#'                     segments = 5, segment.type="consecutive", 
#'                     ncomp=c(16,15))
#'                     
#' # Print results
#' print(pem_NRP)
#' 
#' # Plot all results
#' plot(pem_NRP)
#' }
#' # Reversed order of NIR and Raman (uncomment to run)
#' # paths_RNP <- matrix(c(2,1,2,3,1,3), ncol = 2, byrow = TRUE)
#' # pem_RNP <- path_effects(paths_RNP, emulsions, validation="CV",
#' #                   segments = 5, segment.type="consecutive", ncomp=c(16,15))
#' # print(pem_RNP)
#' @importFrom pls plsr R2
#' @importFrom pracma Rank
#'
#' @export
path_effects <- function(relations, blocks, validation, segments, SO=TRUE, 
                         fits=FALSE, boot = 0, ncomp = NULL, transitive_closure = FALSE, ...){
  # Validate relations structure early
  if(!is.matrix(relations) || ncol(relations) != 2){
    stop("relations must be a 2-column matrix of block indices")
  }
  if(any(!is.finite(relations)) || any(abs(relations - round(relations)) > 0)){
    stop("relations must contain finite integer block indices")
  }
  relations <- matrix(as.integer(relations), ncol = 2)

  # Check that blocks is a data frame or list of matrices
  if (is.data.frame(blocks)) {
    blocks <- as.list(blocks)
  } else if (!is.list(blocks) || !all(sapply(blocks, is.matrix))) {
    stop("blocks must be a data frame or a list of matrices")
  }
  # Check that all blocks have names
  if (is.null(names(blocks)) || any(names(blocks) == "")) {
    stop("All blocks must have names")
  }
  nblock <- length(blocks)
  if(any(relations < 1L | relations > nblock)){
    stop("relations contain block indices outside 1:", nblock)
  }

  # Keep a de-duplicated copy of input relations for ncomp parsing in original block indexing
  relations_input <- unique(relations)

  # Parse ncomp in the original block indexing.
  # It will be re-ordered after topological sorting.
  ncomp_parsed <- NULL
  if(!is.null(ncomp)){
    if(!is.numeric(ncomp)){
      stop("ncomp must be a numeric vector of integers")
    }
    if(any(!is.finite(ncomp))){
      stop("ncomp contains non-finite values")
    }
    if(any(abs(ncomp - round(ncomp)) > 0)){
      stop("ncomp must contain integer values")
    }
    ncomp <- as.integer(ncomp)

    predictor_idx <- sort(unique(relations_input[,1]))
    ncomp_block <- rep(0L, nblock)

    if(is.null(names(ncomp))){
      if(length(ncomp) == length(predictor_idx)){
        ncomp_block[predictor_idx] <- ncomp
      } else if(length(ncomp) == nblock){
        ncomp_block[] <- ncomp
      } else {
        stop("ncomp must have length equal to number of predictor blocks (", length(predictor_idx), ") or total blocks (", nblock, ")")
      }
    } else {
      nm <- names(ncomp)
      idx_from_names <- match(nm, names(blocks))
      if(any(is.na(idx_from_names))){
        suppressWarnings(idx_numeric <- as.integer(nm))
        if(any(is.na(idx_numeric)) || any(idx_numeric < 1 | idx_numeric > nblock)){
          stop("Named ncomp entries must match block names or valid block indices")
        }
        idx_from_names <- idx_numeric
      }
      ncomp_block[idx_from_names] <- ncomp
      if(any(ncomp_block[predictor_idx] == 0L)){
        stop("ncomp must specify non-zero values for all predictor blocks")
      }
    }

    if(any(ncomp_block[predictor_idx] == 0L)){
      stop("ncomp must specify non-zero values for all predictor blocks")
    }
    signs <- sign(ncomp_block[predictor_idx])
    if(length(unique(signs)) > 1){
      stop("ncomp cannot mix positive (limit) and negative (forced) values")
    }
    ncomp_parsed <- ncomp_block
  }

  # Internally reorder blocks by topological order so directed edges such as 2->1
  # are represented in the upper-triangular decomposition loops.
  edges <- unique(relations_input)
  indegree <- integer(nblock)
  children <- vector("list", nblock)
  for(k in seq_len(nrow(edges))){
    from <- edges[k,1]
    to <- edges[k,2]
    indegree[to] <- indegree[to] + 1L
    children[[from]] <- c(children[[from]], to)
  }
  for(v in seq_len(nblock)){
    if(length(children[[v]]) > 0){
      children[[v]] <- sort(unique(children[[v]]))
    }
  }
  queue <- which(indegree == 0L)
  topo <- integer(0)
  while(length(queue) > 0){
    v <- queue[1]
    queue <- queue[-1]
    topo <- c(topo, v)
    if(length(children[[v]]) > 0){
      for(w in children[[v]]){
        indegree[w] <- indegree[w] - 1L
        if(indegree[w] == 0L){
          queue <- c(queue, w)
        }
      }
    }
  }
  if(length(topo) != nblock){
    stop("relations must define an acyclic graph (DAG)")
  }
  remap <- integer(nblock)
  remap[topo] <- seq_len(nblock)
  blocks <- blocks[topo]
  relations <- cbind(remap[edges[,1]], remap[edges[,2]])
  if(!is.null(ncomp_parsed)){
    ncomp <- ncomp_parsed[topo]
  }

  block_nrows <- vapply(blocks, nrow, integer(1))
  if(length(unique(block_nrows)) != 1){
    stop("All blocks must have the same number of rows")
  }
  N <- block_nrows[1]
  
  # Precompute total number of individual variables for blocks that will be processed
  total_ind_cols <- 0
  for(j in 2:nblock){
    include_j <- FALSE
    if(transitive_closure){
      include_j <- TRUE
    } else {
      # Check if any path ends at block j
      include_j <- any(relations[,2]==j)
    }
    if(include_j){
      total_ind_cols <- total_ind_cols + ncol(blocks[[j]])
    }
  }
  
  mat <- matS <- matrix(NA, nrow = (nblock-1)*4, ncol = nblock-1)
  mat_ind <- matS_ind <- matrix(NA, nrow = (nblock-1)*4, ncol = total_ind_cols)
  miss_U <- matrix(FALSE, nrow = (nblock-1), ncol = nblock-1)
  miss_U_ind <- matrix(FALSE, nrow = (nblock-1), ncol = total_ind_cols)
  if(boot > 0){
    mat_boot <- matS_boot <- array(NA, dim = c((nblock-1)*4, nblock-1, boot))
  }
  colStart <- numeric(nblock)  # Pre-allocate colStart vector
  colStart[1] <- 1
  # Pre-compute colStart values and ind_names for all blocks
  ind_names <- character()
  for(j in 2:nblock){
    include_j <- FALSE
    if(transitive_closure){
      include_j <- TRUE
    } else {
      include_j <- any(relations[,2]==j)
    }
    if(include_j){
      colStart[j] <- colStart[j-1] + ncol(blocks[[j]])
      ind_names <- c(ind_names, colnames(blocks[[j]]))
    } else {
      colStart[j] <- colStart[j-1]
    }
  }
  scaleWarn <- scaleWarnInd <- 0
  s_ind <- vector("list", nblock-1)  # Pre-allocate s_ind list
  ncomps <- vector("list", nblock-1)  # Pre-allocate list
  for(i in 1:(nblock-1)){
    for(j in (i+1):nblock){
      if(transitive_closure || any(relations[,1]==i & relations[,2]==j)){
      result <- TotUnCoAd_regression(relations, i, j, blocks, validation = validation, 
                         segments = segments, SO = SO, fits = fits, ncomp = ncomp, ...)
      mat[((i-1)*4+1):((i-1)*4+4),j-1] <- c(result$Tot_EF, result$Un_EF, result$Co_EF, result$Ad_EF)
      if(result$All_EF==0){
        scaleWarn <- 1
        s <- result$Tot_EF
      } else {
        s <- result$All_EF
      }
      if(s==0){
        scaleWarn <- 2
        s <- 1
      }
      matS[((i-1)*4+1):((i-1)*4+4),j-1] <- c(result$Tot_EF, result$Un_EF, result$Co_EF, result$Ad_EF)/s
      if(boot > 0){
        for(b in 1:boot){
          ind <- sample(N,N, replace = TRUE)
          blocks_boot <- lapply(blocks, function(bl) bl[ind, , drop = FALSE])
          result_boot <- TotUnCoAd_regression(relations, i, j, blocks_boot, validation = validation, 
                                  segments = segments, SO = SO, fits = fits, ncomp = ncomp, ...)
          mat_boot[((i-1)*4+1):((i-1)*4+4), j-1, b] <- c(result_boot$Tot_EF, result_boot$Un_EF, result_boot$Co_EF, result_boot$Ad_EF)
          if(result_boot$All_EF==0){
            scaleWarn <- 1
            s <- result_boot$Tot_EF
          } else {
            s <- result_boot$All_EF
          }
          if(s==0){
            scaleWarn <- 2
            s <- 1
          }
          matS_boot[((i-1)*4+1):((i-1)*4+4), j-1, b] <- c(result_boot$Tot_EF, result_boot$Un_EF, result_boot$Co_EF, result_boot$Ad_EF)/s
        }
      }
      # Initialize ncomps[[i]] if not already done
      if(is.null(ncomps[[i]]))
        ncomps[[i]] <- list()
      ncomps[[i]][[j]] <- result$ncomps
      # Initialize per-target scaling vector the first time this target block is computed
      if(is.null(s_ind[[j-1]])){
        if(result$All_EF==0){
          scaleWarnInd <- 1
          s_ind[[j-1]] <- result$Tot_EF_ind
        } else {
          s_ind[[j-1]] <- result$All_EF_ind
        }
        if(any(s_ind[[j-1]]==0)){
          scaleWarnInd <- 2
          s_ind[[j-1]][s_ind[[j-1]]==0] <- 1
        }
      }
      if(colStart[j] > colStart[j-1] && !is.null(s_ind[[j-1]])){  # Only assign if columns were allocated and scaling is available
        mat_ind[((i-1)*4+1):((i-1)*4+4),colStart[j-1]:(colStart[j]-1)] <- rbind(result$Tot_EF_ind, result$Un_EF_ind, result$Co_EF_ind, result$Ad_EF_ind)
        matS_ind[((i-1)*4+1):((i-1)*4+4),colStart[j-1]:(colStart[j]-1)] <- rbind(result$Tot_EF_ind/s_ind[[j-1]], result$Un_EF_ind/s_ind[[j-1]], result$Co_EF_ind/s_ind[[j-1]], result$Ad_EF_ind/s_ind[[j-1]])
      }
      miss_U[i,j-1] <- result$miss_U
      if(colStart[j] > 0){  # Only update if columns were allocated for block j
        miss_U_ind[i,colStart[j-1]:(colStart[j]-1)] <- result$miss_U
      }
      } else {
        miss_U[i,j-1] <- TRUE
        # Only update miss_U_ind if columns were allocated for block j
        if(colStart[j] > 0){
          miss_U_ind[i,colStart[j-1]:(colStart[j]-1)] <- TRUE
        }
      }
    }
  }
  # Add block names to ncomps entries that were computed
  for(i in seq_along(ncomps)){
    if(!is.null(ncomps[[i]]) && length(ncomps[[i]]) > 0){
      # Get j indices that have data (indices are j values since we use ncomps[[i]][[j]])
      idx_with_data <- which(!sapply(ncomps[[i]], is.null))
      # Create names vector with proper names for populated entries
      all_names <- rep("", length(ncomps[[i]]))
      all_names[idx_with_data] <- names(blocks)[idx_with_data]
      names(ncomps[[i]]) <- all_names
    }
  }
  names(ncomps) <- names(blocks)[-nblock]
  
  colnames(mat_ind) <- colnames(matS_ind) <- ind_names
  colnames(mat) <- colnames(matS) <- names(blocks)[-1]
  rownames(mat) <- rownames(matS) <- rownames(mat_ind) <- rownames(matS_ind) <- paste(rep(c("Tot","Un ","Co ","Ad "),nblock-1), rep(names(blocks)[-nblock],each=4), sep=" ")
  class(mat) <- c("path_effects", class(mat))
  attr(mat, "miss_U") <- miss_U
  attr(mat, "scaled") <- matS
  attr(mat, "scaleWarn") <- scaleWarn
  attr(mat, "individual") <- mat_ind
  attr(mat, "individualScaled") <- matS_ind
  attr(mat, "scaleWarnInd") <- scaleWarnInd
  attr(mat, "miss_U_ind") <- miss_U_ind
  attr(mat, "colStart") <- colStart
  attr(mat, "ncomps") <- ncomps
  attr(mat, "transitive_closure") <- transitive_closure
  if(boot > 0){
    attr(mat, "boot") <- mat_boot
    attr(mat, "bootScaled") <- matS_boot
  }
  mat
}

#' Print Path Effects
#'
#' Prints a summary of path effects analysis with optional bootstrap standard errors.
#'
#' @param x An object of class `path_effects`.
#' @param from Optional source block name or index. If specified with `to`, prints only that block pair.
#' @param to Optional target block name or index. If specified with `from`, prints only that block pair.
#' @param nsmall Number of decimal places to display. Default is 2.
#' @param digits Number of significant digits. Default is 1.
#' @param scaled Logical. If `TRUE`, display scaled contributions. Default is `FALSE`.
#' @param individual Logical. If `TRUE`, display individual-level effects. Default is `FALSE`.
#' @param boot Logical. If `TRUE` and bootstrap samples are available, display standard errors. 
#'   Default is `TRUE`.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns `x`.
#'
#' @seealso [path_effects()], [plot.path_effects()]
#'
#' @export
print.path_effects <- function(x, from=NULL, to=NULL, nsmall = 2, digits = 1, scaled = FALSE, individual=FALSE, boot=TRUE, ...){
  x0 <- x
  rn <- rownames(x)
  miss_U <- attr(x, "miss_U")
  print_transitive <- isTRUE(attr(x, "transitive_closure"))
  
  # Handle from/to selection for single block-pair print
  if(xor(is.null(from), is.null(to))){
    stop("Both 'from' and 'to' must be provided together")
  }
  print_single <- !is.null(from) && !is.null(to)
  if(print_single){
    nblock <- ncol(x0) + 1
    source_names <- rownames(x0)[seq(1, by=4, length.out=nblock-1)]
    source_names <- substring(source_names, 5)
    target_names <- colnames(x0)
    if(is.character(from)){
      from_idx <- which(source_names == from)
      if(length(from_idx) == 0) stop("Block '", from, "' not found")
    } else {
      from_idx <- from
    }
    if(is.character(to)){
      to_idx <- which(target_names == to)
      if(length(to_idx) == 0) stop("Block '", to, "' not found")
    } else {
      to_idx <- to
    }
    if(from_idx < 1 || from_idx > nblock-1) stop("Source block index out of range")
    if(to_idx < 1 || to_idx > nblock-1) stop("Target block index out of range")
  }

  if(boot && !is.null(attr(x, "boot")) && !scaled){
    boot <- TRUE
    if(scaled){
      xb <- attr(x, "bootScaled")
    } else {
      xb <- attr(x, "boot")
    }
  } else {
    boot <- FALSE
  }
  if(scaled){
    # Scale by the All_EF part
    if(individual){
      scaleWarn <- attr(x, "scaleWarnInd")
      miss_U <- attr(x, "miss_U_ind")
      x <- attr(x, "individualScaled")
      if(scaleWarn == 1){
        warning("The individual combined contribution is 0, so scaling by the total contribution")
      } else if(scaleWarn == 2){
        warning("The individual total contribution and combined contributions are 0, so scaling by 1")
      }
    } else{
      scaleWarn <- attr(x, "scaleWarn")
      x <- attr(x, "scaled")
      if(scaleWarn == 1){
        warning("The combined contribution is 0, so scaling by the total contribution")
      } else if(scaleWarn == 2){
        warning("The total contribution and combined contributions are 0, so scaling by 1")
      }
    }
  } else {
    if(individual){
      miss_U <- attr(x, "miss_U_ind")
      x <- attr(x, "individual")
    }
  }

  if(print_single){
    row_idx <- ((from_idx-1)*4+1):((from_idx-1)*4+4)
    if(individual){
      colStart <- attr(x0, "colStart")
      col_idx <- colStart[to_idx]:(colStart[to_idx+1]-1)
      x <- x[row_idx, col_idx, drop=FALSE]
      miss_U <- miss_U[from_idx, col_idx, drop=FALSE]
      if(boot){
        # Bootstrap summaries are not available for individual outputs
        boot <- FALSE
      }
    } else {
      x <- x[row_idx, to_idx, drop=FALSE]
      miss_U <- matrix(miss_U[from_idx, to_idx], nrow = 1, ncol = 1)
      if(boot){
        xb <- xb[row_idx, to_idx, , drop=FALSE]
      }
    }
  }
  
  xc <- gsub("NaN","  x  ", gsub("NA","  ", format(x[]*100, digits=digits, nsmall=nsmall, scientific=FALSE)))
  for(i in 1:(nrow(x)/4)){
    for(j in 1:ncol(x)){
      if(miss_U[i,j] && j>=i && !print_transitive){
        xc[(i-1)*4+1,j] <- "  x  "
        xc[(i-1)*4+2,j] <- "  x  "
        xc[(i-1)*4+3,j] <- "  x  "
        xc[(i-1)*4+4,j] <- "  x  "
      }
    }
  }
  if(boot){
    xb <- apply(xb,1:2,sd)
    xbc <- gsub("NaN","  x  ", gsub("NA","  ", format(xb[]*100, digits=digits, nsmall=nsmall, scientific=FALSE)))
    xbc <- matrix(paste0("(",xbc, ")"), nrow=nrow(xbc))
    rownames(xbc) <- rownames(xc)
    colnames(xbc) <- paste0("std(",colnames(xc),")")
    for(i in 1:(nrow(x)/4)){
      for(j in 1:ncol(x)){
        if(miss_U[i,j] && j>=i && !print_transitive){
          xbc[(i-1)*4+1,j] <- ""
          xbc[(i-1)*4+2,j] <- ""
          xbc[(i-1)*4+3,j] <- ""
          xbc[(i-1)*4+4,j] <- ""
        }
        if(j<i){
          xbc[(i-1)*4+1,j] <- ""
          xbc[(i-1)*4+2,j] <- ""
          xbc[(i-1)*4+3,j] <- ""
          xbc[(i-1)*4+4,j] <- ""
        }
      }
    }
    for(i in 1:nrow(x)){
      for(j in 1:ncol(x)){
        if(grepl("x",xc[i,j])){
          xbc[i,j] <- ""
        }
      }
    }
    
    xc <- cbind(xc,xbc)[,c(matrix(1:(ncol(x)*2),2, byrow=TRUE))]
  }
  print(xc, quote=FALSE)
}

#' @keywords internal
#' @noRd
effectplot <- function(Tot_EF, Un_EF, Co_EF, Ad_EF, 
                       ylim = c(-20,100), ylab = "", xlab = "", space=0.2, spec=FALSE, ...){
  if(length(spec)>0){
    plot(spec,Tot_EF, type="l", ylim=ylim, ylab=ylab, xlab=xlab, lwd=2, ...)
    lines(spec,Co_EF, col="black", lwd=1, lty=2)
    lines(spec,Co_EF+Un_EF, col="gray", lwd=1)
    lines(spec,Co_EF+Un_EF+Ad_EF, col="gray", lwd=1, lty=2)
  } else {
    barplot(Ad_EF, offset=Un_EF+Co_EF, col="gray", density=20, border="black", ylim=ylim, 
            ylab=ylab, xlab=xlab, space=space, ...)#, names.arg=c("A","B","C","D"))
    barplot(Un_EF, offset=Co_EF, col="lightgray", border="black", add=TRUE, space=space, ...)
    barplot(Co_EF, col=NA, border="black", density=10, add=TRUE, space=space, ...)
    barplot(100+numeric(length(Un_EF)), border="black", col=NA, add=TRUE, space=space, ...)
    abline(h=0, lwd=2)
    for(i in 1:length(Un_EF)){
      lines(c(space+(i-1)*(1+space)-0.03,(1+space)+(i-1)*(1+space)+0.03),c(Tot_EF[i],Tot_EF[i]), lwd=2)
    }
  }
}

#' @keywords internal
#' @noRd
plot.TotUnCoAd_regression <- function(x, ylab = "", mar = c(0.5,2.8,1.8,0.3), scaled=FALSE, ...){
  nvar <- length(x$Tot_EF_ind)
  if(nvar < 4){
    r <- 1
    c <- nvar
  } else {
    c <- ceiling(sqrt(nvar))
    r <- ceiling(nvar/c)
  }
  par.old <- par(mfrow=c(r, c), mar = mar, mgp=c(1.8,0.5,0))
  for(i in 1:nvar){
    if(scaled){
      scaleWarn <- 0
      if(x$All_EF_ind[i]==0){
        scaleWarn <- 1
        s <- x$Tot_EF_ind[i]
      } else {
        s <- x$All_EF_ind[i]
      }
      if(s==0){
        scaleWarn <- 2
        s <- 1
      }    
      if(scaleWarn == 1){
        warning("The combined contribution is 0, so scaling by the total contribution")
      } else if(scaleWarn == 2){
        warning("The total contribution and combined contributions are 0, so scaling by 1")
      }
    } else {
      s <- 1
    }
    
    Tot_EF <- x$Tot_EF_ind[i]*100/s
    Un_EF <- x$Un_EF_ind[i]*100/s
    Co_EF <- x$Co_EF_ind[i]*100/s
    Ad_EF <- x$Ad_EF_ind[i]*100/s
    names(Tot_EF) <- names(Un_EF) <- names(Co_EF) <- names(Ad_EF) <- NULL
    if(!is.finite(x$Co_EF_ind[i])){
      Co_EF <- rep(0, length(Tot_EF))
      x$Co_EF_ind <- Co_EF
    }
    if(!is.finite(x$Ad_EF_ind[i])){
      Ad_EF <- rep(0, length(Tot_EF))
    }
    effectplot(Tot_EF, Un_EF, Co_EF, Ad_EF, 
               ylim = c(min(min(x$Co_EF_ind),0),100), 
               ylab = ylab, xlab="")
    axis(3, 0.7, x$vars[i], lwd.ticks = 0, mgp=c(1.8,0.5,0))
  }
  par(par.old)
}

#' Plot Path Effects
#'
#' Visualisation of multiblock path effects decomposition.
#'
#' @param x An object of class `path_effects`.
#' @param from Optional source block name or index. If specified with `to`, plots only that block pair.
#' @param to Optional target block name or index. If specified with `from`, plots only that block pair.
#' @param scaled Logical. If `TRUE`, display scaled contributions. Default is `FALSE`.
#' @param individual Logical. If `TRUE`, display individual-level effects. Default is `FALSE`.
#' @param spectra Optional numeric vector for spectral data visualization. Default is empty.
#' @param mar Margins for the plot. Default is `c(0.5,2.8,1.8,0.3)`.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return Invisibly returns `NULL`.
#'
#' @seealso [path_effects()], [print.path_effects()]
#'
#' @export
plot.path_effects <- function(x, from=NULL, to=NULL, scaled=FALSE, individual=FALSE, spectra=numeric(0), mar = c(0.5,2.8,1.8,0.3), ...){
  nblock <- ncol(x)+1
  miss_U <- attr(x, "miss_U")
  miss_U_ind <- attr(x, "miss_U_ind")
  colStart <- attr(x, "colStart")
  vars <- colnames(attr(x,"individual"))
  plot_transitive <- isTRUE(attr(x, "transitive_closure"))
  
  # Handle from/to parameters for single block pair plotting
  plot_single <- !is.null(from) && !is.null(to)
  if(plot_single){
    # Convert block names to indices if needed
    source_names <- rownames(x)[seq(1, by=4, length.out=nblock-1)]
    source_names <- substring(source_names, 5)  # Remove "T  " prefix
    target_names <- colnames(x)
    if(is.character(from)){
      i <- which(source_names == from)
      if(length(i)==0) stop("Block '", from, "' not found")
    } else {
      i <- from
    }
    if(is.character(to)){
      j <- which(target_names == to)
      if(length(j)==0) stop("Block '", to, "' not found")
    } else {
      j <- to
    }
    if(i < 1 || i > nblock-1) stop("Source block index out of range")
    if(j < 1 || j > nblock-1) stop("Target block index out of range")
  }
  
  if(scaled){
    # Scale by the All_EF part
    if(individual){
      scaleWarn <- attr(x, "scaleWarnInd")
      x <- attr(x, "individualScaled")
      if(scaleWarn == 1){
        warning("The individual combined contribution is 0, so scaling by the total contribution")
      } else if(scaleWarn == 2){
        warning("The individual total contribution and combined contributions are 0, so scaling by 1")
      }
    } else {
      scaleWarn <- attr(x, "scaleWarn")
      x <- attr(x, "scaled")
      if(scaleWarn == 1){
        warning("The combined contribution is 0, so scaling by the total contribution")
      } else if(scaleWarn == 2){
        warning("The total contribution and combined contributions are 0, so scaling by 1")
      }
    }
  } else {
    if(individual)
      x <- attr(x, "individual")
  }
  x[!is.finite(x)] <- 0
  
  if(plot_single){
    # Single block pair plot
    par.old <- par(mar = mar, mgp=c(1.8,0.5,0))
    if(miss_U[i,j] && !plot_transitive){
      warning(
        paste0(
          "Path '", substring(rownames(x)[(i-1)*4+1], 5), " -> ", colnames(x)[j],
          "' is not available with transitive_closure = FALSE. ",
          "Re-run path_effects(..., transitive_closure = TRUE) to plot this path."
        )
      )
      plot(0,0, type="n", axes=FALSE, xlab="", ylab="")
    } else {
      if(individual){
        Tot_EF <- x[(i-1)*4+1,colStart[j]:(colStart[j+1]-1), drop=FALSE]*100
        Un_EF <- x[(i-1)*4+2,colStart[j]:(colStart[j+1]-1), drop=FALSE]*100
        Co_EF <- x[(i-1)*4+3,colStart[j]:(colStart[j+1]-1), drop=FALSE]*100
        Ad_EF <- x[(i-1)*4+4,colStart[j]:(colStart[j+1]-1), drop=FALSE]*100
      } else {
        Tot_EF <- x[(i-1)*4+1,j]*100
        Un_EF <- x[(i-1)*4+2,j]*100
        Co_EF <- x[(i-1)*4+3,j]*100
        Ad_EF <- x[(i-1)*4+4,j]*100
      }
      if(length(spectra) > 0 && (j+1) %in% spectra){
        spec <- as.numeric(vars[colStart[j]:(colStart[j+1]-1)])
      } else {
        spec <- numeric(0)
      }
      effectplot(Tot_EF, Un_EF, Co_EF, Ad_EF, 
                 ylim = c(min(min(x),0),100), ylab=substring(rownames(x)[(i-1)*4+1], 5), 
                 xlab=colnames(x)[j], spec=spec, ...)
    }
    par(par.old)
  } else {
    # Full grid plot (original behavior)
    par.old <- par(mfrow=c(nblock-1, nblock-1), mar = mar, mgp=c(1.8,0.5,0))
    for(i in 1:(nblock-1)){
      for(j in 1:(nblock-1)){
        if(j < i){
          plot(0,0, type="n", axes=FALSE, xlab="", sub = "", ylab="")
        } else {
          if(miss_U[i,j] && !plot_transitive){
            plot(0,0, type="n", axes=FALSE, xlab="", sub = "", ylab="")
          } else {
            if(individual){
              Tot_EF <- x[(i-1)*4+1,colStart[j]:(colStart[j+1]-1), drop=FALSE]*100
              Un_EF <- x[(i-1)*4+2,colStart[j]:(colStart[j+1]-1), drop=FALSE]*100
              Co_EF <- x[(i-1)*4+3,colStart[j]:(colStart[j+1]-1), drop=FALSE]*100
              Ad_EF <- x[(i-1)*4+4,colStart[j]:(colStart[j+1]-1), drop=FALSE]*100
            } else {
              Tot_EF <- x[(i-1)*4+1,j]*100
              Un_EF <- x[(i-1)*4+2,j]*100
              Co_EF <- x[(i-1)*4+3,j]*100
              Ad_EF <- x[(i-1)*4+4,j]*100
            }
            if(length(spectra) > 0 && (j+1) %in% spectra){
              spec <- as.numeric(vars[colStart[j]:(colStart[j+1]-1)])
            } else {
              spec <- numeric(0)
            }
            effectplot(Tot_EF, Un_EF, Co_EF, Ad_EF, 
                       ylim = c(min(min(x),0),100), ylab=substring(rownames(x)[(i-1)*4+1], 5), xlab="", spec=spec, ...)
            if(!individual){
              axis(3, 0.7, paste0(colnames(x)[j],ifelse(miss_U[i,j],"*","")), lwd.ticks = 0, mgp=c(1.8,0.5,0))
            }
          }
        }
      }
    }
    par(par.old)
  }
}
