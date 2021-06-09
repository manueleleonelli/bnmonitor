#' Frobenius norm
#'
#' \code{Fro} returns the Frobenius norm between a Bayesian network and its update after parameter variation.
#'
#' The details depend on the class the method \code{Fro} is applied to.
#'
#' @seealso \code{\link{KL.GBN}}, \code{\link{KL.CI}}, \code{\link{Fro.CI}}, \code{\link{Fro.GBN}}, \code{\link{Jeffreys.GBN}}, \code{\link{Jeffreys.CI}}
#' @param x object of class \code{GBN} or \code{CI}.
#' @param ... parameters specific to the class used.
#' @export
#'
#'@return A dataframe whose columns depend of the class of the object.
#'
Fro <- function (x, ...) {
  UseMethod("Fro", x)
}

#' Frobenius norm for \code{GBN}
#'
#' \code{Fro.GBN} returns the Frobenius norm between between an object of class \code{GBN}  and its update after a standard parameter variation.
#'
#' Computation of the Frobenius norm between a Bayesian network and the additively perturbed Bayesian network, where the perturbation is either to the mean vector or to the covariance matrix. The Frobenius norm is not computed for perturbations of the mean since it is always equal to zero.
#'
#'@param x object of class \code{GBN}.
#'@param entry a vector of length 2 indicating the entry of the covariance matrix to vary.
#'@param delta numeric vector, including the variation parameters that act additively.
#'@param log boolean value. If \code{TRUE}, the logarithm of the Frobenius norm is returned. Set by default to \code{TRUE}.
#'@param ... additional arguments for compatibility.
#'
#'@return A dataframe including in the first column the variations performed and in the second column the corresponding Frobenius norm.
#'
#'@examples Fro(synthetic_gbn,c(3,3),seq(-1,1,0.1))
#'
#' @seealso \code{\link{KL.GBN}}, \code{\link{KL.CI}}, \code{\link{Fro.CI}}, \code{\link{Jeffreys.GBN}}, \code{\link{Jeffreys.CI}}
#'
#' @export
#'@importFrom ggplot2 ggplot
#'@importFrom ggplot2 geom_line
#'@importFrom ggplot2 geom_point
#'@importFrom ggplot2 labs
#'@importFrom ggplot2 aes
#'
Fro.GBN <- function(x,entry,delta, log = TRUE, ...){
  gbn <- x
  fro <- numeric(length(delta))
    D <- matrix(0,length(gbn$mean),length(gbn$mean))
    for(i in 1:length(fro)){
      D[entry[1],entry[2]]<- delta[i]
      D[entry[2],entry[1]]<- delta[i]
      if(is.psd(round(gbn$covariance+D,2))){
        fro[i] <- sum(diag(t(D)%*%D))
      }
      else{fro[i]<-NA}
    }
  fro <- data.frame(Variation = delta,Frobenius=fro)
  if(log == TRUE){fro[,-1] <- log(fro[,-1])}
  Variation <- fro$Variation
  Frobenius <- fro$Frobenius
    if(nrow(fro)==1){
      plot <- suppressWarnings(ggplot(data = fro, mapping = aes(x = Variation, y = Frobenius)) + geom_point( na.rm = T) + labs(x = "delta",  y = "Frobenius", title = "Frobenius norm") + theme_minimal())
    }else{
      plot <- suppressWarnings(ggplot(data = fro, mapping = aes(x = Variation, y = Frobenius)) + geom_line( na.rm = T) + labs(x = "delta",  y = "Frobenius", title = "Frobenius norm") + theme_minimal())
    }
  out <- list(Frobenius = fro, plot = plot)
  attr(out,'class') <- 'fro'
  return(out)
}

#' Frobenius norm for \code{CI}
#'
#' \code{Fro.CI} returns the Frobenius norm between an object of class \code{CI}  and its update after a model-preserving parameter variation.
#'
#' Computation of the Frobenius norm between a Bayesian network and its updated version after a model-preserving variation.
#'
#'@param x object of class \code{CI}.
#'@param type character string. Type of model-preserving co-variation: either \code{"total"}, \code{"partial"}, \code{row}, \code{column} or \code{all}. If \code{all} the Frobenius norm is computed for every type of co-variation matrix.
#'@param entry a vector of length 2 indicating the entry of the covariance matrix to vary.
#'@param delta numeric vector with positive elements, including the variation parameters that act multiplicatively.
#'@param log boolean value. If \code{TRUE}, the logarithm of the Frobenius norm is returned. Set by default to \code{TRUE}.
#'@param ... additional arguments for compatibility.
#'
#'
#'@return A dataframe including in the first column the variations performed, and in the following columns the corresponding Frobenius norms for the chosen model-preserving co-variations.
#'
#'@seealso \code{\link{KL.GBN}}, \code{\link{KL.CI}}, \code{\link{Fro.GBN}}, \code{\link{Jeffreys.GBN}}, \code{\link{Jeffreys.CI}}
#'@references C. Görgen & M. Leonelli (2020), Model-preserving sensitivity analysis for families of Gaussian distributions.  Journal of Machine Learning Research, 21: 1-32.
#'
#'@examples Fro(synthetic_ci,"total",c(1,1),seq(0.9,1.1,0.01))
#'@examples Fro(synthetic_ci,"partial",c(1,4),seq(0.9,1.1,0.01))
#'@examples Fro(synthetic_ci,"column",c(1,2),seq(0.9,1.1,0.01))
#'@examples Fro(synthetic_ci,"row",c(3,2),seq(0.9,1.1,0.01))
#'
#'@export
#'

Fro.CI <- function(x, type, entry, delta, log = TRUE, ...){
  ci <- x
  fro <- numeric(length(delta))
  if(type == "total"){
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov <- total_covar_matrix(ci,entry,delta[i])
      if(is.psd(round(Cov*Delta*ci$covariance,2))){
        fro[i] <- sum(diag(t(ci$covariance-Cov*Delta*ci$covariance)%*%(ci$covariance-Cov*Delta*ci$covariance)))
      }
      else{fro[i]<- NA}
    }
  }
  if(type == "partial"){
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov <- partial_covar_matrix(ci,entry,delta[i])
      if(is.psd(round(Cov*Delta*ci$covariance,2))){
        fro[i] <- sum(diag(t(ci$covariance-Cov*Delta*ci$covariance)%*%(ci$covariance-Cov*Delta*ci$covariance)))
      }
      else{fro[i]<- NA}
    }
  }
  if(type == "row"){
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov <- row_covar_matrix(ci, entry, delta[i])
      if(is.psd(round(Cov*Delta*ci$covariance,2))){
        fro[i] <- sum(diag(t(ci$covariance-Cov*Delta*ci$covariance)%*%(ci$covariance-Cov*Delta*ci$covariance)))      }
      else{fro[i]<- NA}
    }
  }
  if(type == "column"){
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov <- col_covar_matrix(ci,entry,delta[i])
      if(is.psd(round(Cov*Delta*ci$covariance,2))){
        fro[i] <- sum(diag(t(ci$covariance-Cov*Delta*ci$covariance)%*%(ci$covariance-Cov*Delta*ci$covariance)))      }
      else{fro[i]<- NA}
    }
    }
  if(type == "all"){
    fro <- matrix(0,length(delta),4)
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov_col <- col_covar_matrix(ci,entry,delta[i])
      Cov_row <- row_covar_matrix(ci,entry,delta[i])
      Cov_par <- partial_covar_matrix(ci,entry,delta[i])
      Cov_tot <- total_covar_matrix(ci,entry,delta[i])
      if(is.psd(round(Cov_tot*Delta*ci$covariance,2))){
        fro[i,1] <- sum(diag(t(ci$covariance-Cov_tot*Delta*ci$covariance)%*%(ci$covariance-Cov_tot*Delta*ci$covariance)))      }
      else{fro[i,1]<- NA}
      if(is.psd(round(Cov_par*Delta*ci$covariance,2))){
        fro[i,2] <- sum(diag(t(ci$covariance-Cov_par*Delta*ci$covariance)%*%(ci$covariance-Cov_par*Delta*ci$covariance)))      }
      else{fro[i,2]<- NA}
      if(is.psd(round(Cov_row*Delta*ci$covariance,2))){
        fro[i,3] <- sum(diag(t(ci$covariance-Cov_row*Delta*ci$covariance)%*%(ci$covariance-Cov_row*Delta*ci$covariance)))      }
      else{fro[i,3]<- NA}
      if(is.psd(round(Cov_col*Delta*ci$covariance,2))){
        fro[i,4] <- sum(diag(t(ci$covariance-Cov_col*Delta*ci$covariance)%*%(ci$covariance-Cov_col*Delta*ci$covariance)))      }
      else{fro[i,4]<- NA}
    }
  }
  if(type == "all")
    {
    Fro_data <- data.frame(Variation = delta, Total = fro[,1], Partial = fro[,2], Row_based = fro[,3], Column_based = fro[,4])
    }
  else{
    Fro_data <- data.frame(Variation = delta, Frobenius= fro)
    }
  if(log == TRUE){Fro_data[,-1] <- log(Fro_data[,-1])}
    if(type == "all"){
    ci <- gather(Fro_data, "scheme", "value", - Variation)
    scheme <- ci$scheme
    value <- ci$value
    Variation <- ci$Variation
      if(nrow(Fro_data) == 1){
        plot <-  ggplot(data = ci, mapping = aes(x = Variation, y = value)) + geom_point(aes(color = scheme)) + labs( x = "delta", y = "Frobenius", title = "Frobenius Norm") + theme_minimal()
      } else{
        plot <-  ggplot(data = ci, mapping = aes(x = Variation, y = value)) + geom_line(aes(color = scheme)) + labs(x = "delta",  y = "Frobenius", title = "Frobenius Norm") + theme_minimal()
      }
    } else{
      Variation <- Fro_data$Variation
      Frobenius <- Fro_data$Frobenius
      if(nrow(Fro_data) == 1){
        plot <- suppressWarnings(ggplot(data = Fro_data, mapping = aes(x = Variation, y = Frobenius)) + geom_point( na.rm = T) + labs(x = "delta", y = "Frobenius", title = "Frobenius norm") + theme_minimal())
      }else{
        plot <- suppressWarnings(ggplot(data = Fro_data, mapping = aes(x = Variation, y = Frobenius)) + geom_line(na.rm = T ) + labs(x = "delta",  y = "Frobenius", title = "Frobenius norm") + theme_minimal())
      }
    }
  out <- list(Frobenius = Fro_data, plot = plot)
  attr(out,'class') <- 'fro'
  return(out)

}




is.psd <- function (x, tol = 1e-08)
{
  if (!is.square.matrix(x))
    stop("argument x is not a square matrix")
  if (!is.symmetric.matrix(x))
    stop("argument x is not a symmetric matrix")
if (!is.numeric(x))
  stop("argument x is not a numeric matrix")
eigenvalues <- eigen(x, only.values = TRUE)$values
n <- nrow(x)
for (i in 1:n) {
  if (abs(eigenvalues[i]) < tol) {
    eigenvalues[i] <- 0
  }
}
if (any(eigenvalues < 0)) {
  return(FALSE)
}
return(TRUE)
}

is.square.matrix <- function (x)
{
  if (!is.matrix(x))
    stop("argument x is not a matrix")
  return(nrow(x) == ncol(x))
}

is.symmetric.matrix <- function (x)
{
  if (!is.matrix(x)) {
    stop("argument x is not a matrix")
  }
  if (!is.numeric(x)) {
    stop("argument x is not a numeric matrix")
  }
  if (!is.square.matrix(x))
    stop("argument x is not a square numeric matrix")
  return(sum(x == t(x)) == (nrow(x)^2))
}

