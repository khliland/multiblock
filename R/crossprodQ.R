crossprodQ <- function(A,B){
  if(missing(B)){
    if(storage.mode(A) == "double"){
      return(crossprodnum(A))
    }
    if(storage.mode(A) == "integer"){
      return(crossprodint(A))
    }
    if(storage.mode(A) == "logical"){
      storage.mode(A) <- "integer"
      return(crossprodint(A))
    }
  } else {
    a <- storage.mode(A)
    b <- storage.mode(B)
    if(a==b){
      if(a=="double"){
        return(crossprodnumnum(A,B))
      }
      if(a=="integer"){
        return(crossprodintint(A,B))
      }
      if(a=="logical"){
        storage.mode(A) <- "integer"
        storage.mode(B) <- "integer"
        return(crossprodintint(A,B))
      }
    } else {
      if(any(c(a,b)%in%"double")){
        storage.mode(A) <- "double"
        storage.mode(B) <- "double"
        return(crossprodnumnum(A,B))
      } else {
        storage.mode(A) <- "integer"
        storage.mode(B) <- "integer"
        return(crossprodintint(A,B))
      }
    }
  }
}

tcrossprodQ <- function(A,B){
  if(missing(B)){
    if(storage.mode(A) == "double"){
      return(tcrossprodnum(A))
    }
    if(storage.mode(A) == "integer"){
      return(tcrossprodint(A))
    }
    if(storage.mode(A) == "logical"){
      storage.mode(A) <- "integer"
      return(tcrossprodint(A))
    }
  } else {
    a <- storage.mode(A)
    b <- storage.mode(B)
    if(a==b){
      if(a=="double"){
        return(tcrossprodnumnum(A,B))
      }
      if(a=="integer"){
        return(tcrossprodintint(A,B))
      }
      if(a=="logical"){
        storage.mode(A) <- "integer"
        storage.mode(B) <- "integer"
        return(tcrossprodintint(A,B))
      }
    } else {
      if(any(c(a,b)%in%"double")){
        storage.mode(A) <- "double"
        storage.mode(B) <- "double"
        return(tcrossprodnumnum(A,B))
      } else {
        storage.mode(A) <- "integer"
        storage.mode(B) <- "integer"
        return(tcrossprodintint(A,B))
      }
    }
  }
}
