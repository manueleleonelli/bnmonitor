#' Standard variation of the mean vector
#'
#' Computation of an updated \code{GBN} object after a variation of the mean vector.
#'
#' Need to write this description, write examples and explain the return
#'
#' @family gbn sensitivity
#'
#'@param gbn object of class \code{GBN}.
#'@param entry an index specifying the entry of the mean vector to vary.
#'@param delta additive variation coefficient for the entry of the mean vector given in \code{entry}.
#'
#'@references Gómez-Villegas, M. A., Maín, P., & Susi, R. (2007). Sensitivity analysis in Gaussian Bayesian networks using a divergence measure. Communications in Statistics—Theory and Methods, 36(3), 523-539.
#'@references Gómez-Villegas, M. A., Main, P., & Susi, R. (2013). The effect of block parameter perturbations in Gaussian Bayesian networks: Sensitivity and robustness. Information Sciences, 222, 439-458.
#'@export


mean_var <- function(gbn,entry,delta){
  if (length(entry)!= 1){stop("entry must be a positive integer")}
  if (entry>length(gbn$mean) || entry <= 0 || entry%%1!=0){stop("entry must be a positive integer")}
  gbn$mean[entry] <- gbn$mean[entry] + delta
  return(gbn)
}

#' Standard variation of the covariance matrix
#'
#' Computation of an updated \code{GBN} object after a variation of the covariance matrix.
#'
#' Need to write this description, write examples and explain the return
#'
#' @family gbn sensitivity
#'
#'@param gbn object of class \code{GBN}.
#'@param entry a vector of length 2 specifying the entry of the covariance matrix to vary.
#'@param delta additive variation coefficient for the entry of the mean vector given in \code{entry}.
#'
#'@references Gómez-Villegas, M. A., Maín, P., & Susi, R. (2007). Sensitivity analysis in Gaussian Bayesian networks using a divergence measure. Communications in Statistics—Theory and Methods, 36(3), 523-539.
#'@references Gómez-Villegas, M. A., Main, P., & Susi, R. (2013). The effect of block parameter perturbations in Gaussian Bayesian networks: Sensitivity and robustness. Information Sciences, 222, 439-458.
#' @importFrom matrixcalc is.positive.semi.definite
#'@export


covariance_var <- function(gbn, entry, delta){
  if (delta <= 0){stop("delta must be strictly positive")}
  if (length(entry)!= 2){stop("entry must be a vector of length 2")}
  if (entry[1]>nrow(gbn$covariance) || entry[1] <= 0 || entry[1]%%1!=0){stop("rows of entry not appropriate")}
  if (entry[2]>ncol(gbn$covariance) || entry[2] <= 0 || entry[1]%%1!=0){stop("columns of entry not appropriate")}
  var_matrix <- variation_mat(gbn,entry,delta)
  gbn$covariance <- var_matrix + gbn$covariance
  if(is.positive.semi.definite(gbn$covariance)){return(gbn)}
  else{
    gbn$warning <- "The covariance is not positive semidefinite"
    attr(gbn,'class') <- 'npsd.gbn'
    return(gbn)}
}
