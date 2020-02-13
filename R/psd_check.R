#' Check for positive semidefiniteness after a perturbation
#'
#' \code{psd_check} returns a boolean to determine if the covariance matrix after a perturbation is positive semidefinite.
#'
#' The details depend on the class the method \code{KL} is applied to.
#'
#' @param x object of class \code{GBN} or \code{CI}.
#' @param ... parameters specific to the class used.
#'
#'
#' @seealso \code{\link{psd_check.GBN}}, \code{\link{psd_check.CI}}
#'
#' @export
#'

psd_check <- function (x, ...) {
  UseMethod("psd_check", x)
}

#' Check for positive semidefiniteness after a perturbation for \code{GBN}
#'
#' \code{psd_check} returns a boolean to determine if the covariance matrix after a perturbation is positive semidefinite.
#'
#' Let \eqn{\Sigma} be the covariance matrix of a Gaussian Bayesian network and let \eqn{D} be a perturbation matrix acting additively. The perturbed covariance matrix \eqn{\Sigma+D} is positive semidefinite if
#' \deqn{\rho(D)\leq \lambda_{\min}(\Sigma)}
#' where \eqn{\lambda_{\min}} is the smallest eigenvalue end \eqn{\rho} is the spectral radius.
#'
#'@seealso \code{\link{psd_check.CI}}
#'
#'@param x object of class \code{GBN}.
#'@param entry a vector of length 2 indicating the entry of the covariance matrix to vary.
#'@param delta numeric vector, including the variation parameters that act additively.
#'@param ... additional arguments for compatibility.
#'
#'@return A dataframe including the variations performed and the check for positive semidefinitess.
#'
#'@references Goergen, C., & Leonelli, M. (2018). Model-preserving sensitivity analysis for families of Gaussian distributions. arXiv preprint arXiv:1809.10794.
#'
#'@examples psd_check(synthetic_gbn,c(2,4),-3)
#'@examples psd_check(synthetic_gbn,c(2,3),seq(-1,1,0.1))
#'
#'@export
#'
#'

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


#' Check for positive semidefiniteness after a perturbation for \code{CI}
#'
#' \code{psd_check} returns a boolean to determine if the covariance matrix after a perturbation is positive semidefinite.
#'
#' Let \eqn{\Sigma} be the covariance matrix of a Gaussian Bayesian network, \eqn{\Delta} be a perturbation matrix acting multiplicatively and \eqn{\tilde\Delta} a model-preserving covariation. The perturbed covariance matrix \eqn{\tilde\Delta\circ\Delta\circ\Sigma} is positive semidefinite if
#' \deqn{\rho(\tilde\Delta\circ\Delta\circ\Sigma-\Sigma)\leq \lambda_{\min}(\Sigma)}
#' where \eqn{\lambda_{\min}} is the smallest eigenvalue, \eqn{\rho} is the spectral radius and \eqn{\circ} is the Schur or elementwise product.
#'
#'@seealso \code{\link{psd_check.GBN}}
#'
#'@param x object of class \code{CI}.
#'@param type character string. Type of model-preserving covariation: either \code{"total"}, \code{"partial"}, \code{row}, \code{column} or \code{all}. If \code{all} the check is performed for every type of covariation matrix.
#'@param entry a vector of length 2 indicating the entry of the covariance matrix to vary.
#'@param delta numeric vector, including the variation parameters that act multiplicatively.
#'@param ... additional arguments for compatibility.
#'
#'@return A dataframe including the variations performed and the check for positive semidefinitess.
#'
#'@references Goergen, C., & Leonelli, M. (2018). Model-preserving sensitivity analysis for families of Gaussian distributions. arXiv preprint arXiv:1809.10794.
#'
#'@examples psd_check(synthetic_ci,"partial",c(2,4),0.95)
#'@examples psd_check(synthetic_ci,"all",c(2,3),seq(0.9,1.1,0.01))
#'@export
#'
#'

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
