#' Model-Preserving covariation
#'
#' Model-Preserving covariation for objects of class \code{CI}.
#'
#' Need to write this description, write examples and explain the return
#'
#'@seealso \code{\link{covariance_var}}, \code{\link{covariation_matrix}}
#'
#'@param ci object of class \code{CI}.
#'@param type character string. Type of model-preserving covariation: either \code{"total"}, \code{"partial"}, \code{row} or \code{column}.
#'@param entry a vector of length two specifying the entry of the covariance matrix to vary.
#'@param delta multiplicative variation coefficient for the entry of the covariance matrix given in \code{entry}.
#'
#'@references Goergen, C., & Leonelli, M. (2018). Model-preserving sensitivity analysis for families of Gaussian distributions. arXiv preprint arXiv:1809.10794.
#'@importFrom matrixcalc is.positive.semi.definite
#'@export


model_pres_cov <- function(ci, type, entry, delta){
  if (delta <= 0){stop("delta must be strictly positive")}
  if (length(entry)!= 2){stop("entry must be a vector of length 2")}
  if (type != "total" & type != "partial" & type != "row" & type != "column"){stop("wrong covariation scheme")}
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
    if(is.positive.semi.definite(ci$covariance)){return(ci)}
    else{
      ci$warning <- "The covariance is not positive semidefinite"
      attr(ci,'class') <- 'npsd.ci'
      return(ci)
    }
  }
  else{
    ci$covariance <- var_matrix*ci$covariance
    if(is.positive.semi.definite(ci$covariance)){return(ci)}
    else{
      ci$warning <- "The covariance is not positive semidefinite"
      attr(ci,'class') <- 'npsd.ci'
      return(ci)
    }
  }
}

