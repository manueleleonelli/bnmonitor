#' Model-Preserving co-variation
#'
#' Model-preserving co-variation for objects of class \code{CI}.
#'
#' Let the original Bayesian network have a Normal distribution \eqn{\mathcal{N}(\mu,\Sigma)} and let \code{entry} be equal to \eqn{(i,j)}. For a multiplicative variation of the covariance matrix by an amount \eqn{\delta}, a variation matrix \eqn{\Delta} is constructed as
#' \deqn{\Delta_{k,l}=\left\{
#' \begin{array}{ll}
#' \delta & \mbox{if } k=i, l=j\\
#' \delta & \mbox{if } l=i, k=j \\
#' 0 & \mbox{otherwise}
#' \end{array}
#' \right.}
#' A co-variation matrix \eqn{\tilde\Delta} is then constructed and the resulting distribution after the variation is \eqn{\mathcal{N}(\mu,\tilde\Delta\circ\Delta\circ\Sigma)}, assuming \eqn{\tilde\Delta\circ\Delta\circ\Sigma} is positive semi-definite and where \eqn{\circ} denotes the Schur (or element-wise) product. The matrix \eqn{\tilde\Delta} is so constructed to ensure that all conditional independence in the original Bayesian networks are retained after the parameter variation.


#'@seealso \code{\link{covariance_var}}, \code{\link{covariation_matrix}}
#'
#'
#'@examples model_pres_cov(synthetic_ci,"partial",c(1,3),1.1)
#'@examples model_pres_cov(synthetic_ci,"partial",c(1,3),0.9)
#'@examples model_pres_cov(synthetic_ci,"total",c(1,2),0.5)
#'@examples model_pres_cov(synthetic_ci,"row",c(1,3),0.98)
#'@examples model_pres_cov(synthetic_ci,"column",c(1,3),0.98)
#'
#'
#'@param ci object of class \code{CI}.
#'@param type character string. Type of model-preserving co-variation: either \code{"total"}, \code{"partial"}, \code{row} or \code{column}.
#'@param entry a vector of length two specifying the entry of the covariance matrix to vary.
#'@param delta multiplicative variation coefficient for the entry of the covariance matrix given in \code{entry}.
#'
#'@return If the resulting covariance is positive semi-definite, \code{model_pres_cov} returns an object of class \code{CI} with an updated covariance matrix. Otherwise it returns an object of class \code{npsd.ci}, which has the same components of \code{CI} but also has a warning entry specifying that the covariance matrix is not positive semi-definite.
#'
#'@references C. Görgen & M. Leonelli (2020), Model-preserving sensitivity analysis for families of Gaussian distributions.  Journal of Machine Learning Research, 21: 1-32.
#'@export


model_pres_cov <- function(ci, type, entry, delta){
  if (delta <= 0){stop("delta must be strictly positive")}
  if (length(entry)!= 2){stop("entry must be a vector of length 2")}
  if (type != "total" & type != "partial" & type != "row" & type != "column"){stop("wrong co-variation scheme")}
  var_matrix <- variation_mat(ci,entry,delta)
  test <- rep(F,length(ci$cond_ind))
  for(i in 1:length(ci$cond_ind)){
    if(ci$order[entry[1]] %in% unique(c(ci$cond_ind[[i]]$A,ci$cond_ind[[i]]$C)) & ci$order[entry[2]] %in% unique(c(ci$cond_ind[[i]]$B,ci$cond_ind[[i]]$C))){test[i] <- T}
    if(ci$order[entry[2]] %in% unique(c(ci$cond_ind[[i]]$A,ci$cond_ind[[i]]$C)) & ci$order[entry[1]] %in% unique(c(ci$cond_ind[[i]]$B,ci$cond_ind[[i]]$C))){test[i] <- T}
  }
  if(sum(test)>0){
    if(type == "total"){cov_matrix <- total_covar_matrix(ci,entry,delta)}
    if(type == "partial"){cov_matrix <- partial_covar_matrix(ci,entry,delta)}
    if(type == "row"){cov_matrix <- row_covar_matrix(ci,entry, delta)}
    if(type == "column"){cov_matrix <- col_covar_matrix(ci,entry,delta)}
    ci$covariance <- cov_matrix*var_matrix*ci$covariance
    if(is.psd(round(ci$covariance,5))){return(ci)}
    else{
      ci$warning <- "The covariance is not positive semidefinite"
      attr(ci,'class') <- 'npsd.ci'
      return(ci)
    }
  }
  else{
    ci$covariance <- var_matrix*ci$covariance
    if(is.psd(round(ci$covariance,5))){return(ci)}
    else{
      ci$warning <- "The covariance is not positive semidefinite"
      attr(ci,'class') <- 'npsd.ci'
      return(ci)
    }
  }
}

