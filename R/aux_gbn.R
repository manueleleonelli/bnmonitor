
simple_ci <- function(ci){
  b.temp <- c()
  for(i in 1:(length(ci[[4]]))){
    if(identical(ci[[4]][[i]]$C,character(0)) || identical(ci[[4]][[i]]$B,character(0))){b.temp <- c(b.temp,F)} else{b.temp <- c(b.temp,T)}
  }
  temp <- vector(mode = "list", length = sum(b.temp))
  count <- 1
  for(i in 1:length(ci[[4]])){
    if (b.temp[i] == T){
      temp[[count]] <- ci$cond_ind[[i]]
      count <- count+1
    }
  }
  ci$cond_ind <- temp
  return(ci)
}

ci_submatrix <- function(ci){
  simpleci <- simple_ci(ci)
  A <- c()
  B <- c()
  C <- c()
  for(i in 1:length(simpleci[[4]])){
    A <- c(A,simpleci[[4]][[i]]$A)
    B <- c(B,simpleci[[4]][[i]]$B)
    C <- c(C,simpleci[[4]][[i]]$C)
  }
  return(list(A=unique(A),B=unique(B),C=unique(C)))
}



variation_mat <- function(x,entry,delta){
  mat <- matrix(1,nrow = nrow(x$covariance),ncol= ncol(x$covariance))
  mat[entry[1],entry[2]] <- delta
  mat[entry[2],entry[1]] <- delta
  return(mat)
}

