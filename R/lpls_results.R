#' @name lpls_results
#' @title Result functions for L-PLS objects (\code{lpls})
#' @aliases plot.lpls predict.lpls lplsCV
#' 
#' @description Correlation loading plot, prediction and cross-validation for L-PLS
#' models with class \code{\link{lpls}}.
#'
#' @param x \code{lpls} object
#' @param comps \code{integer} vector of components.
#' @param doplot \code{logical} indicating if plotting should be performed.
#' @param level \code{integer} vector of length 3 for selecting plot symbol. 1=dots. 2=dimnames.
#' @param arrow \code{integer} vector of length 3 indicating arrows (1) or not (0).
#' @param xlim \code{numeric} x limits.
#' @param ylim \code{numeric} y limits.
#' @param samplecol \code{character} for sample colours.
#' @param pathcol \code{character} for third colour.
#' @param varcol \code{character} for variable colours.
#' @param varsize \code{numeric} size of symbols for variables.
#' @param sampleindex \code{integer} for selecting samples.
#' @param pathindex \code{integer} for selecting in third direction.
#' @param varindex \code{integer} for selecting variables.
#' @param object \code{lpls} object.
#' @param X1new \code{matrix} of new X1 samples.
#' @param X2new \code{matrix} of new X2 samples.
#' @param X3new \code{matrix} of new X3 samples.
#' @param exo.direction \code{character} selecting "X2" or "X3" prediction.
#' @param segments1 \code{list} of sample segments.
#' @param segments2 \code{list} of variable segments.
#' @param trace \code{logical} indicating if verbose mode should be selected.
#' @param ... Not implemented.
#'
#' @return Nothing is return for plotting (\code{plot.lpls}), predicted values are returned for predictions (\code{predict.lpls})
#' and cross-validation metrics are returned for for cross-validation (\code{lplsCV}).
#'
#' @examples
#' # Simulate data set
#' sim <- lplsData(I = 30, N = 20, J = 5, K = 6, ncomp = 2)
#' X1  <- sim$X1; X2 <- sim$X2; X3 <- sim$X3
#' 
#' # exo-L-PLS:
#' lp.exo  <- lpls(X1,X2,X3, ncomp = 2)
#' # Predict X1
#' pred.exo.X2 <- predict(lp.exo, X1new = X1, exo.direction = "X2")
#' # Predict X3
#' pred.exo.X2 <- predict(lp.exo, X1new = X1, exo.direction = "X3")
#' 
#' # endo-L-PLS:
#' lp.endo <- lpls(X1,X2,X3, ncomp = 2, type = "endo")
#' # Predict X1 from X2 and X3 (in this case fitted values):
#' pred.endo.X1 <- predict(lp.endo, X2new = X2, X3new = X3)
#' 
#' # LOO cross-validation horizontally
#' lp.cv1 <- lplsCV(lp.exo, segments1 = as.list(1:dim(X1)[1]))
#' 
#' # LOO cross-validation vertically
#' lp.cv2 <- lplsCV(lp.exo, segments2 = as.list(1:dim(X1)[2]))
#' 
#' # Three-fold CV, horizontal
#' lp.cv3 <- lplsCV(lp.exo, segments1 = as.list(1:10, 11:20, 21:30))
#' 
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' @export
plot.lpls <- function(x, comps=c(1,2), doplot=c(TRUE,TRUE,TRUE), level=c(2,2,2),
                      arrow=c(1,0,1), xlim=c(-1,1), ylim=c(-1,1), samplecol=4, pathcol=2, varcol="grey70",
                      varsize=1, sampleindex=1:dim(x$corloadings$R22)[1], pathindex=1:dim(x$corloadings$R3)[1],
                      varindex=1:dim(x$corloadings$R21)[1], ...){
  
  plottype <- c("p","n")
  
  plot(xlim,ylim,"n",
       xlab=paste("Comp", comps[1]),
       ylab=paste("Comp", comps[2]),
       main=paste("Correlation loading plot from ",x$options$type,"-L-PLS analysis",sep=""))
  
  ellipse(c(0,0), matrix(c(1,0,0,1),2,2), radius=1, lty=1, lwd=1, col=1, center.pch=FALSE)
  abline(h=0, lty=3)
  abline(v=0, lty=3)   
  
  if(doplot[1]){
    points(x$corloadings$R21[varindex,comps[1]], x$corloadings$R21[varindex,comps[2]],type=plottype[level[2]],pch=20,col=varcol,cex=varsize)
    points(x$corloadings$R22[sampleindex,comps[1]], x$corloadings$R22[sampleindex,comps[2]],type=plottype[level[2]],pch=20,col=samplecol,cex=2)
    if(level[2]==2){
      text(x$corloadings$R21[varindex,comps[1]], x$corloadings$R21[varindex,comps[2]],labels=colnames(x$data$X2),cex=varsize,col=varcol)
      text(x$corloadings$R22[sampleindex,comps[1]], x$corloadings$R22[sampleindex,comps[2]],labels=rownames(x$data$X2),cex=0.7,col=samplecol)
    }
    if(arrow[2]==1){
      arrows(0,0,x$corloadings$R21[varindex,comps[1]], x$corloadings$R21[varindex,comps[2]],lwd=2,col="grey70",length = 0.1)
      arrows(0,0,x$corloadings$R22[sampleindex,comps[1]], x$corloadings$R22[sampleindex,comps[2]],lwd=2,col=4,length = 0.1)    
    }
  }
  
  if(doplot[2]){
    points(x$corloadings$R1[,comps[1]],x$corloadings$R1[,comps[2]],type=plottype[level[1]],pch=20,col=3,cex=2)
    if(level[1]==2){
      text(x$corloadings$R1[,comps[1]],x$corloadings$R1[,comps[2]],labels=colnames(x$data$X1),cex=0.7,col=3)    
    }
    if(arrow[1]==1){arrows(0,0,x$corloadings$R1[,comps[1]],x$corloadings$R1[,comps[2]],col=3,length = 0.1)}
  }
  
  if(doplot[3]){
    points(x$corloadings$R3[pathindex,comps[1]],x$corloadings$R3[pathindex,comps[2]],type=plottype[level[3]],pch=20,col=2,cex=2)
    if(level[3]==2){
      text(x$corloadings$R3[pathindex,comps[1]],x$corloadings$R3[pathindex,comps[2]],labels=colnames(x$data$X3),cex=0.7,col=pathcol)    
    }
    if(arrow[3]==1){arrows(0,0,x$corloadings$R3[pathindex,comps[1]],x$corloadings$R3[pathindex,comps[2]],col=2,length = 0.1)}
  }
}

#' @rdname lpls_results 
#' @export
predict.lpls <- function(object, X1new = NULL, X2new = NULL, X3new = NULL, exo.direction = c("X2", "X3"), ...){

  # Rename inputs from Smilde, Liland and Næs 2021 to Sæbø et al.
  X1newa <- X1new
  X1new  <- X2new
  X2new  <- X1newa
  if(!is.null(X3new))
    X3new  <- t(X3new)
  
  # Prediction for endo-LPLS
  if(object$options$type=="endo"){
    
    if(any(is.na(X1new)) | any(is.na(X3new))) stop("Prediction requires complete predictor data.\n In case of cross-validation, use 'impute=TRUE' in model fit\n")
    k <- dim(X1new)[2]
    ntest <-dim(X1new)[1]
    l <- dim(X3new)[2]
    ptest <- dim(X3new)[1]
    
    # Centering and scaling of new observations
    if(object$options$scale[2]){
      X1new <- scale(X1new,object$means$mX1,attr(object$data$X1,"scaled:scale"))
    } else {
      X1new <- scale(X1new,object$means$mX1, scale=F)
    }
    
    if(object$options$scale[3]){  
      X3new <- scale(X3new,object$means$mX3, attr(object$data$X3,"scaled:scale"))
    } else {
      X3new <- scale(X3new,object$means$mX3, scale=F)
    }
    
    pred <- matrix(1,nrow=ntest,ncol=ptest)*object$means$grandmX2 + X1new%*%object$coefficients$C%*%t(X3new)
    if(!is.null(attr(object$data$X2,"scaled:scale"))) cat("Warning: Only grand mean adjusted and column-scaled X2 has been predicted")
    
    #Not implemented yet. Correction for scaling and double centering.
    #    if(is.null(attr(object$data$X2,"scaled:scale"))){
    #      pred <- matrix(1,nrow=ntest,ncol=ptest)*object$means$grandmX2 + X1new%*%object$coefficients$C%*%t(X3new)
    #    }else{
    #      pred <- matrix(1,nrow=ntest,ncol=ptest)*object$means$grandmX2 + X1new%*%object$coefficients$C%*%t(X3new)*(matrix(1,ntest,1)%*%attr(object$data$X2,"scaled:scale"))
    #    }
    
    
    res <- list(pred=pred)
    
  } else if(object$options$type=="exo"){
    n <- dim(object$data$X2)[1]
    p <- dim(object$data$X2)[2]    
    if(any(is.na(X2new))) stop("Prediction requires complete predictor data.\n In case of cross-validation, use 'impute=TRUE' in model fit\n")
    X2new <- as.matrix(X2new)
    ntest <- dim(X2new)[1]
    ptest <- dim(X2new)[2]
    rowm1 <- apply(X2new,1,mean)
    colm1 <- apply(X2new,2,mean)
    
    if(exo.direction=="X2"){   #Prediction of X1 (Sæbø name)
      
      # Double centering of new obs
      if(!object$options$doublecenter){
        X2new <- X2new - matrix(object$means$grandmX2,nrow=ntest,ncol=p)
      } else {
        X2new <- X2new-
          t(matrix(rep(1,p),ncol=1)%*%rowm1) - 
          matrix(rep(1,ntest),ncol=1)%*%object$means$colmX2 + 
          matrix(object$means$grandmX2,nrow=ntest,ncol=p)    
      }

      if(is.null(attr(object$data$X1,"scaled:scale"))){
        pred <- matrix(1,ntest,1)%*%object$means$mX1 + X2new%*%object$coefficients$B1
      } else {
        pred <- matrix(1,ntest,1)%*%object$means$mX1 + X2new%*%object$coefficients$B1*(matrix(1,ntest,1)%*%attr(object$data$X1,"scaled:scale"))
      }
      res <- list(pred=pred)
    }
    
    if(exo.direction=="X3"){   #Prediction of X3
      
      # Centering of new obs
      if(!object$options$doublecenter){
        X2new <- X2new - matrix(object$means$grandmX2,nrow=n,ncol=ptest)
      } else {
        X2new <- X2new-
          t(matrix(rep(1,ptest),ncol=1)%*%object$means$rowmX2) - 
          matrix(rep(1,n),ncol=1)%*%colm1 + 
          matrix(object$means$grandmX2,nrow=n,ncol=ptest)    
      }
      
      if(is.null(attr(object$data$X3,"scaled:scale"))){
        pred <- matrix(1,ptest,1)%*%object$means$mX3 + t(X2new)%*%object$coefficients$B3
      } else {
        pred <- matrix(1,ptest,1)%*%object$means$mX3 + t(X2new)%*%object$coefficients$B3*(matrix(1,ptest,1)%*%attr(object$data$X3,"scaled:scale"))
      }                
      
      res <- list(pred = pred)
    }
    
  } else if(object$options$type=="exo_ort"){
    cat("No prediction method is implemented for method exo_ort\n")
    res <- NULL
  }
  res
}

#' @rdname lpls_results 
#' @export
lplsCV <- function(object, segments1 = NULL, segments2 = NULL, trace = TRUE){
  X1 <- object$data$X1
  X2 <- object$data$X2
  X3 <- object$data$X3
  n <- dim(object$data$X2)[1]
  if(is.null(segments1) & is.null(segments2)) segments1 <- as.list(1:n)
  
  if(object$options$type=="exo_ort") stop("No prediction method for type='exo_ort'. Use type='exo'.\n")
  if(object$options$type=="exo"){
    if(!is.null(segments1)&!is.null(segments2)) stop("CV can only be run in one direction for type='exo'")
    if(!is.null(segments1)){
      nsegs <- length(segments1)
      
      rowpred <- dim(X1)[1]
      colpred <- dim(X1)[2]
      pred <- array(0,dim=c(rowpred, colpred,object$ncomp))
      for(i in 1:nsegs){
        testX1 <- X1[segments1[[i]],,drop=F]
        testX2 <- X2[segments1[[i]],,drop=F]
        for(j in 1:object$ncomp){
          trainfit <- lpls(X2, X1, t(X3), ncomp=j, 
                           doublecenter=object$options$doublecenter,
                           scaledata=object$options$scale,
                           type=object$options$type,
                           subsetX2=-segments1[[i]])
          pred[segments1[[i]],,j] <- predict(trainfit,X1new=testX2, exo.direction="X2")$pred 
        }
        if(trace)cat(paste("Segment",i,"of",nsegs,"completed\n"))
      }
      dimnames(pred) <- list(dimnames(X1)[[1]],dimnames(X1)[[2]],paste("Comp",1:object$ncomp,sep=""))
      rmsep<- apply(pred,3,function(x){sqrt(mean((x-X1)^2))})
    } else {
      nsegs <- length(segments2)
      
      rowpred <- dim(X3)[1]
      colpred <- dim(X3)[2]
      pred <- array(0,dim=c(rowpred, colpred,object$ncomp))
      
      for(i in 1:nsegs){
        testX3 <- X3[segments2[[i]],,drop=F]
        testX2 <- X2[,segments2[[i]],drop=F]
        
        for(j in 1:object$ncomp){
          trainfit <- lpls(X2, X1, t(X3), ncomp=j, 
                           doublecenter=object$options$doublecenter,
                           scaledata=object$options$scale,
                           type=object$options$type,
                           subsetX3=-segments2[[i]])
          pred[segments2[[i]],,j] <- predict(trainfit, X1new=testX2, exo.direction="X3")$pred
        }
        if(trace) cat(paste("Segment",i,"of",nsegs,"completed\n"))
      }
      
      dimnames(pred) <- list(dimnames(X3)[[1]],dimnames(X3)[[2]],paste("Comp",1:object$ncomp,sep=""))
      rmsep<- apply(pred,3,function(x){sqrt(mean((x-X3)^2))})
      
    }
    
  } else if(object$options$type=="endo"){
    if(!is.null(segments1)&!is.null(segments2)) stop("CV can only be run in one direction for type='exo'")
    if(!is.null(segments1)){
      nsegs <- length(segments1)
      
      rowpred <- dim(X2)[1]
      colpred <- dim(X2)[2]
      pred <- array(0,dim=c(rowpred, colpred,object$ncomp))
      for(i in 1:nsegs){
        testX1 <- X1[segments1[[i]],,drop=F]
        testX2 <- X2[segments1[[i]],,drop=F]
        for(j in 1:object$ncomp){
          trainfit <- lpls(X2, X1, t(X3), ncomp=j, 
                           doublecenter=object$options$doublecenter,
                           scaledata=object$options$scale,
                           type=object$options$type,
                           subsetX2=-segments1[[i]])
          pred[segments1[[i]],,j] <- predict(trainfit, X2new=testX1, X3new=t(X3))$pred 
        }
        if(trace) cat(paste("Segment",i,"of",nsegs,"completed\n"))
      }
      dimnames(pred) <- list(dimnames(X2)[[1]],dimnames(X2)[[2]],paste("Comp",1:object$ncomp,sep=""))
      rmsep<- apply(pred,3,function(x){sqrt(mean((x-X2)^2))})
    } else {
      nsegs <- length(segments2)
      
      rowpred <- dim(X2)[1]
      colpred <- dim(X2)[2]
      pred <- array(0,dim=c(rowpred, colpred,object$ncomp))
      for(i in 1:nsegs){
        testX3 <- X3[segments2[[i]],,drop=F]
        testX2 <- X2[,segments2[[i]],drop=F]
        for(j in 1:object$ncomp){
          trainfit <- lpls(X2, X1, t(X3), ncomp=j, 
                           doublecenter=object$options$doublecenter,
                           scaledata=object$options$scale,
                           type=object$options$type,
                           subsetX3=-segments2[[i]])
          pred[,segments2[[i]],j] <- predict(trainfit,X1new=X1,X3new=testX3)$pred
        }
        if(trace) cat(paste("Segment",i,"of",nsegs,"completed\n"))
      }
      dimnames(pred) <- list(dimnames(X2)[[1]],dimnames(X2)[[2]],paste("Comp",1:object$ncomp,sep=""))
      rmsep<- apply(pred,3,function(x){sqrt(mean((x-X2)^2))})
    }
  }
  names(rmsep)<-paste("Comp",1:object$ncomp,sep="")
  return(list(rmsep=rmsep, pred=pred))
}
