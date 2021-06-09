#' Standard variation of the mean vector
#'
#' Computation of an updated \code{GBN} object after a variation of the mean vector.
#'
#' Let the original Bayesian network have a Normal distribution \eqn{\mathcal{N}(\mu,\Sigma)} and let \code{entry} be equal to \eqn{i}. Let \eqn{\mu_i} be the i-th entry of \eqn{\mu}. For a variation of the mean  by an amount \eqn{\delta} the resulting distribution is \eqn{\mathcal{N}(\mu',\Sigma)}, where \eqn{\mu'} is equal to \eqn{\mu} except for the i-th entry which is equal to \eqn{\mu+\delta}.

#'@return An object of class \code{GBN} with an updated mean vector.
#'
#'@examples mean_var(synthetic_gbn,2,3)
#'
#'@param gbn object of class \code{GBN}.
#'@param entry an index specifying the entry of the mean vector to vary.
#'@param delta additive variation coefficient for the entry of the mean vector given in \code{entry}.
#'
#'@seealso \code{\link{covariance_var}}
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
#' Let the original Bayesian network have a Normal distribution \eqn{\mathcal{N}(\mu,\Sigma)} and let \code{entry} be equal to \eqn{(i,j)}. For a variation of the covariance matrix by an amount \eqn{\delta}, a variation matrix \eqn{D} is constructed as
#' \deqn{D_{k,l}=\left\{
#' \begin{array}{ll}
#' \delta & \mbox{if } k=i, l=j\\
#' \delta & \mbox{if } l=i, k=j \\
#' 0 & \mbox{otherwise}
#' \end{array}
#' \right.}
#' Then the resulting distribution after the variation is \eqn{\mathcal{N}(\mu,\Sigma +D)}, assuming \eqn{\Sigma+ D} is positive semi-definite.
#'
#'@return If the resulting covariance is positive semi-definite, \code{covariance_var} returns an object of class \code{GBN} with an updated covariance matrix. Otherwise it returns an object of class \code{npsd.gbn}, which has the same components of \code{GBN} but also has a warning entry specifying that the covariance matrix is not positive semi-definite.
#'
#'@param gbn object of class \code{GBN}.
#'@param entry a vector of length 2 specifying the entry of the covariance matrix to vary.
#'@param delta additive variation coefficient for the entry of the co-variation matrix given in \code{entry}.
#'
#'@seealso \code{\link{mean_var}}, \code{\link{model_pres_cov}}
#'
#'@references Gómez-Villegas, M. A., Maín, P., & Susi, R. (2007). Sensitivity analysis in Gaussian Bayesian networks using a divergence measure. Communications in Statistics—Theory and Methods, 36(3), 523-539.
#'@references Gómez-Villegas, M. A., Main, P., & Susi, R. (2013). The effect of block parameter perturbations in Gaussian Bayesian networks: Sensitivity and robustness. Information Sciences, 222, 439-458.
#'
#'@examples covariance_var(synthetic_gbn,c(1,1),3)
#'@examples covariance_var(synthetic_gbn,c(1,2),-0.4)
#'
#'@export


covariance_var <- function(gbn, entry, delta){
  if (length(entry)!= 2){stop("entry must be a vector of length 2")}
  if (entry[1]>nrow(gbn$covariance) || entry[1] <= 0 || entry[1]%%1!=0){stop("rows of entry not appropriate")}
  if (entry[2]>ncol(gbn$covariance) || entry[2] <= 0 || entry[1]%%1!=0){stop("columns of entry not appropriate")}
  if(entry[1]==entry[2] & gbn$covariance[entry[1],entry[2]] + delta < 0){stop("a variance cannot be negative")}
  D <- matrix(0,length(gbn$mean),length(gbn$mean))
    D[entry[1],entry[2]]<- delta
    D[entry[2],entry[1]]<- delta
  gbn$covariance <- D + gbn$covariance
  if(is.psd(round(gbn$covariance,5))){return(gbn)}
  else{
    gbn$warning <- "The covariance is not positive semidefinite"
    attr(gbn,'class') <- 'npsd.gbn'
    return(gbn)}
}
