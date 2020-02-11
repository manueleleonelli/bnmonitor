#' Frobenius norm
#'
#' \code{Fro} returns the Frobenius norm between a Bayesian network and its update after parameter variation.
#'
#' The details depend on the class the method \code{Fro} is applied to.
#'
#' @seealso \code{\link{KL.GBN}}, \code{\link{KL.CI}}, \code{\link{Fro.CI}}, \code{\link{Fro.GBN}}
#' @param x object of class \code{GBN} or \code{CI}.
#' @param ... parameters specific to the class used.
#' @export
#'
#'
Fro <- function (x, ...) {
  UseMethod("Fro", x)
}

#' Frobenius norm for \code{GBN}
#'
#' \code{Fro.GBN} returns the Frobenius norm between a Gaussian Bayesian network and its update after an additive parameter variation.
#'
#'
#'@param x object of class \code{GBN}.
#'@param entry a vector of length 2 indicating the entry of the covariance matrix to vary.
#'@param delta variation parameter that acts additively.
#'@param ... additional parameters to be added to the plot.
#'
#' @seealso \code{\link{KL.GBN}}, \code{\link{KL.CI}}, \code{\link{Fro.CI}}
#'
#' @export
#' @importFrom matrixcalc is.positive.semi.definite
Fro.GBN <- function(x,entry,delta, ...){
  gbn <- x
  fro <- numeric(length(delta))
    D <- matrix(0,length(gbn$mean),length(gbn$mean))
    for(i in 1:length(fro)){
      D[entry[1],entry[2]]<- delta[i]
      D[entry[2],entry[1]]<- delta[i]
      if(is.positive.semi.definite(gbn$covariance+D)){
        fro[i] <- sum(diag(t(D)%*%D))
      }
      else{fro[i]<-NA}
    }
  return(data.frame(Variation = delta,Frobenius=fro))
}

#' Frobenius norm for \code{CI}
#'
#' \code{Fro.CI} returns the Frobenius norm between a Gaussian Bayesian network and its update after a model-preserving parameter variation.
#'
#'@param x object of class \code{CI}.
#'@param type character string: either \code{mean} or \code{covariance} for variations of the mean vector and covariance matrix respectively.
#'@param entry a vector of length 2 indicating the entry of the covariance matrix to vary.
#'@param delta variation parameter that acts multiplicatively.
#'@param ... additional parameters to be added to the plot.
#'
#'@seealso \code{\link{KL.GBN}}, \code{\link{KL.CI}}, \code{\link{Fro.GBN}}
#'@references Goergen, C., & Leonelli, M. (2018). Model-preserving sensitivity analysis for families of Gaussian distributions. arXiv preprint arXiv:1809.10794.
#'
#'@importFrom matrixcalc is.positive.semi.definite
#'@export
#'

Fro.CI <- function(x, type, entry, delta, ...){
  ci <- x
  fro <- numeric(length(delta))
  if(type == "total"){
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov <- total_covar_matrix(ci,entry,delta[i])
      if(is.positive.semi.definite(round(Cov*Delta*ci$covariance,2))){
        fro[i] <- 0.5*(log(det(ci$covariance)/det(Cov*Delta*ci$covariance))-nrow(ci$covariance)+sum(diag(solve(ci$covariance)%*%(Cov*Delta*ci$covariance))))
      }
      else{fro[i]<- NA}
    }
    return(data.frame(Variation = delta, Frobenius = fro))
  }
  if(type == "partial"){
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov <- partial_covar_matrix(ci,entry,delta[i])
      if(is.positive.semi.definite(round(Cov*Delta*ci$covariance,2))){
        fro[i] <- 0.5*(log(det(ci$covariance)/det(Cov*Delta*ci$covariance))-nrow(ci$covariance)+sum(diag(solve(ci$covariance)%*%(Cov*Delta*ci$covariance))))
      }
      else{fro[i]<- NA}
    }
    return(data.frame(Variation = delta, Frobenius = fro))
  }
  if(type == "row"){
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov <- row_covar_matrix(ci, entry, delta[i])
      if(is.positive.semi.definite(round(Cov*Delta*ci$covariance,2))){
        fro[i] <- 0.5*(log(det(ci$covariance)/det(Cov*Delta*ci$covariance))-nrow(ci$covariance)+sum(diag(solve(ci$covariance)%*%(Cov*Delta*ci$covariance))))
      }
      else{fro[i]<- NA}
    }
    return(data.frame(Variation = delta, Frobenius = fro))
  }
  if(type == "column"){
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov <- col_covar_matrix(ci,entry,delta[i])
      if(is.positive.semi.definite(round(Cov*Delta*ci$covariance,2))){
        fro[i] <- 0.5*(log(det(ci$covariance)/det(Cov*Delta*ci$covariance))-nrow(ci$covariance)+sum(diag(solve(ci$covariance)%*%(Cov*Delta*ci$covariance))))
      }
      else{fro[i]<- NA}
    }
    return(data.frame(Variation = delta, Frobenius = fro))
  }
  if(type == "all"){
    fro <- matrix(0,length(delta),4)
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov_col <- col_covar_matrix(ci,entry,delta[i])
      Cov_row <- row_covar_matrix(ci,entry,delta[i])
      Cov_par <- partial_covar_matrix(ci,entry,delta[i])
      Cov_tot <- total_covar_matrix(ci,entry,delta[i])
      if(is.positive.semi.definite(round(Cov_tot*Delta*ci$covariance,2))){
        fro[i,1] <- 0.5*(log(det(ci$covariance)/det(Cov_tot*Delta*ci$covariance))-nrow(ci$covariance)+sum(diag(solve(ci$covariance)%*%(Cov_tot*Delta*ci$covariance))))
      }
      else{fro[i,1]<- NA}
      if(is.positive.semi.definite(round(Cov_par*Delta*ci$covariance,2))){
        fro[i,2] <- 0.5*(log(det(ci$covariance)/det(Cov_par*Delta*ci$covariance))-nrow(ci$covariance)+sum(diag(solve(ci$covariance)%*%(Cov_par*Delta*ci$covariance))))
      }
      else{fro[i,2]<- NA}
      if(is.positive.semi.definite(round(Cov_row*Delta*ci$covariance,2))){
        fro[i,3] <- 0.5*(log(det(ci$covariance)/det(Cov_row*Delta*ci$covariance))-nrow(ci$covariance)+sum(diag(solve(ci$covariance)%*%(Cov_row*Delta*ci$covariance))))
      }
      else{fro[i,3]<- NA}
      if(is.positive.semi.definite(round(Cov_col*Delta*ci$covariance,2))){
        fro[i,4] <- 0.5*(log(det(ci$covariance)/det(Cov_col*Delta*ci$covariance))-nrow(ci$covariance)+sum(diag(solve(ci$covariance)%*%(Cov_col*Delta*ci$covariance))))
      }
      else{fro[i,4]<- NA}
    }
    return(data.frame(Variation = delta, Frobenius_total = fro[,1], Frobenius_partial = fro[,2], Frobenius_row = fro[,3], Frobenius_column = fro[,4]))
  }
}
