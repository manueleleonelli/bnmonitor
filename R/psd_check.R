#' Check for positive semi-definiteness after a perturbation
#'
#' \code{psd_check} returns a boolean to determine if the covariance matrix after a perturbation is positive semi-definite.
#'
#' The details depend on the class the method \code{psd_check} is applied to.
#'
#' Let \eqn{\Sigma} be the covariance matrix of a Gaussian Bayesian network and let \eqn{D} be a perturbation matrix acting additively. The perturbed covariance matrix \eqn{\Sigma+D} is positive semi-definite if
#' \deqn{\rho(D)\leq \lambda_{\min}(\Sigma)}
#' where \eqn{\lambda_{\min}} is the smallest eigenvalue end \eqn{\rho} is the spectral radius.
#'
#'
#' @param x object of class \code{GBN} or \code{CI}.
#' @param type character string. Type of model-preserving co-variation: either \code{total}, \code{partial}, \code{row}, \code{column} or \code{all}. If \code{all}, the Frobenius norms are computed for every type of co-variation matrix.
#'@param entry a vector of length 2 indicating the entry of the covariance matrix to vary.
#'@param delta numeric vector, including the variation parameters that act additively.
#'@param ... additional arguments for compatibility.
#'
#'@return A dataframe including the variations performed and the check for positive semi-definiteness.
#'
#'@references C. Görgen & M. Leonelli (2020), Model-preserving sensitivity analysis for families of Gaussian distributions.  Journal of Machine Learning Research, 21: 1-32.
#'
#'@examples psd_check(synthetic_gbn,c(2,4),-3)
#'@examples psd_check(synthetic_gbn,c(2,3),seq(-1,1,0.1))
#'@examples psd_check(synthetic_ci,"partial",c(2,4),0.95)
#'@examples psd_check(synthetic_ci,"all",c(2,3),seq(0.9,1.1,0.01))
#
#'
#' @export
#'

psd_check <- function (x, ...) {
  UseMethod("psd_check", x)
}


#'@describeIn psd_check \code{psd_check} for objects \code{GBN}
#' @export

psd_check.GBN <- function (x,entry,delta,...){
 gbn <- x
 bol <- rep(TRUE,length(delta))
 min_eig <- min(eigen(gbn$covariance)$values)
 for(i in 1:length(delta)){
 D <- matrix(0,length(gbn$mean),length(gbn$mean))
 D[entry[1],entry[2]] <- delta[i]
 D[entry[2],entry[1]] <- delta[i]
 if(max(abs(eigen(D)$values)) <= min_eig){bol[i] <- TRUE} else{bol[i] <- FALSE}
 }
 return(data.frame(Variation = delta, PSD = bol))
}


#'@describeIn psd_check \code{psd_check} for objects \code{CI}
#' @export

psd_check.CI <- function (x,type,entry,delta,...){
  ci <- x
  bol <- rep(TRUE,length(delta))
  min_eig <- min(eigen(ci$covariance)$values)
  if (type == "partial"){
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov <- partial_covar_matrix(ci,entry,delta[i])
      if(max(abs(eigen(Cov*Delta*ci$covariance-ci$covariance)$values)) <= min_eig){bol[i] <- TRUE} else{bol[i] <- FALSE}
    }
  }
  if (type == "row"){
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov <- row_covar_matrix(ci,entry,delta[i])
      if(max(abs(eigen(Cov*Delta*ci$covariance-ci$covariance)$values)) <= min_eig){bol[i] <- TRUE} else{bol[i] <- FALSE}
    }
  }
  if (type == "column"){
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov <- col_covar_matrix(ci,entry,delta[i])
      if(max(abs(eigen(Cov*Delta*ci$covariance-ci$covariance)$values)) <= min_eig){bol[i] <- TRUE} else{bol[i] <- FALSE}
    }
  }
  if (type == "all"){
    bol <- matrix(TRUE,length(delta),4)
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov_par <- partial_covar_matrix(ci,entry,delta[i])
      Cov_col <- col_covar_matrix(ci,entry,delta[i])
      Cov_row <- row_covar_matrix(ci,entry,delta[i])
      if(max(abs(eigen(Cov_par*Delta*ci$covariance-ci$covariance)$values)) <= min_eig){bol[i,2] <- TRUE} else{bol[i,2] <- FALSE}
      if(max(abs(eigen(Cov_row*Delta*ci$covariance-ci$covariance)$values)) <= min_eig){bol[i,3] <- TRUE} else{bol[i,3] <- FALSE}
      if(max(abs(eigen(Cov_col*Delta*ci$covariance-ci$covariance)$values)) <= min_eig){bol[i,4] <- TRUE} else{bol[i,4] <- FALSE}
      }
  }
  if(type == "all") {return(data.frame(Variation = delta, Total = bol[,1], Partial = bol[,2], Row_based = bol[,3], Column_based = bol[,4]))}else{return(data.frame(Variation = delta, psd.check = bol))}
}
