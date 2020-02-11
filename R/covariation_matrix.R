#' Covariation matrices
#'
#' Construction of model-preserving matrices for objects of class \code{CI}.
#'
#' @seealso \code{\link{model_pres_cov}}
#'@param ci object of class \code{CI}.
#'@param entry a vector of length two specifying the entry of the covariance matrix to vary.
#'@param delta multiplicative variation coefficient for the entry of the covariance matrix given in \code{entry}.
#'
#'@references Goergen, C., & Leonelli, M. (2018). Model-preserving sensitivity analysis for families of Gaussian distributions. arXiv preprint arXiv:1809.10794.
#' @name covariation_matrix
NULL

#' @rdname covariation_matrix
#'
#' @export
total_covar_matrix <- function(ci,entry,delta){
covar <- matrix(delta,nrow = nrow(ci$covariance), ncol = ncol(ci$covariance))
covar[entry[1],entry[2]] <- 1
covar[entry[2],entry[1]] <- 1
return(covar)
}

#' @rdname covariation_matrix
#'
#' @export
col_covar_matrix <- function(ci,entry,delta){
  covar <- matrix(1,nrow = nrow(ci$covariance), ncol = ncol(ci$covariance))
  simpleci<- simple_ci(ci)
  submatrix <- ci_submatrix(simpleci)
  rows <- unique(c(submatrix$A,submatrix$C))
  cols <- unique(c(submatrix$B,submatrix$C))
  temp <- matrix(1,nrow = length(rows), ncol = length(cols))
  ind <- unique(c(which(cols == ci$order[entry[1]]),which(cols == ci$order[entry[2]])))
  ind <- ind[!is.na(ind)]
  while(length(ind)>0){
    for(k in 1:length(ind)){
      for (i in 1:length(rows)){
        temp[i,ind[k]] <- delta
        if(rows[i] %in% cols & cols[ind[k]] %in% rows) {temp[which(cols[ind[k]]==rows),which(cols == rows[i])] <- delta}
      }
    }
    check <- c()
    ind <- c()
    for (i in 1: length(cols)){
      check[i] <- sum(temp[,i]== delta)
      if(check[i] > 0 & check[i] < length(rows)){ind <- c(ind,i)}
    }
  }
  for(i in 1:nrow(covar)){
    for(j in 1:ncol(covar)){
      ind <- c(which(ci$order[i]==rows),which(ci$order[j]==cols))
      if(length(ind) == 2){
        if(temp[ind[1],ind[2]] == delta){
          covar[i,j] <- delta; covar[j,i] <- delta
        }
      }
    }
  }
  return(covar)
}

#' @rdname covariation_matrix
#'
#' @export
partial_covar_matrix <- function(ci, entry, delta){
  covar <- matrix(1,nrow = nrow(ci$covariance), ncol = ncol(ci$covariance))
  simpleci<- simple_ci(ci)
  submatrix <- ci_submatrix(simpleci)
  rows <- unique(c(submatrix$A,submatrix$C))
  cols <- unique(c(submatrix$B,submatrix$C))
  for(i in 1:length(ci$order)){
    for(j in 1:length(ci$order)){
      if(ci$order[i] %in% rows & ci$order[j] %in% cols){
        covar[i,j] <- delta
        covar[j,i] <- delta
      }
    }
  }
  covar[entry[1],entry[2]] <- 1
  covar[entry[2],entry[1]] <- 1
  return(covar)
}


#' @rdname covariation_matrix
#'
#' @export
row_covar_matrix <- function(ci, entry, delta){
  covar <- matrix(1,nrow = nrow(ci$covariance), ncol = ncol(ci$covariance))
  simpleci<- simple_ci(ci)
  submatrix <- ci_submatrix(simpleci)
  rows <- unique(c(submatrix$A,submatrix$C))
  cols <- unique(c(submatrix$B,submatrix$C))
  temp <- matrix(1,nrow = length(rows), ncol = length(cols))
  ind <- unique(c(which(rows == ci$order[entry[1]]),which(rows == ci$order[entry[2]])))
  ind <- ind[!is.na(ind)]
  while(length(ind)>0){
    for(k in 1:length(ind)){
      for (i in 1:length(cols)){
        temp[ind[k],i] <- delta
        if(cols[i] %in% rows & rows[ind[k]] %in% cols) {temp[which(rows == cols[i]),which(rows[ind[k]]==cols)] <- delta}
      }
    }
    check <- c()
    ind <- c()
    for (i in 1: length(rows)){
      check[i] <- sum(temp[i,]== delta)
      if(check[i] > 0 & check[i] < length(cols)){ind <- c(ind,i)}
    }
  }
  for(i in 1:nrow(covar)){
    for(j in 1:ncol(covar)){
      ind <- c(which(ci$order[i]==rows),which(ci$order[j]==cols))
      if(length(ind) == 2){
        if(temp[ind[1],ind[2]] == delta){
          covar[i,j] <- delta; covar[j,i] <- delta
        }
      }
    }
  }
  return(covar)
}
