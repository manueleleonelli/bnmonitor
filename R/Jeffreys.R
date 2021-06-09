#' Jeffreys Divergence
#'
#' \code{Jeffreys} returns the Jeffreys divergence between a continuous Bayesian network and its update after parameter variation.
#'
#' The details depend on the class the method \code{Jefrreys} is applied to.
#'
#' @param x object of class \code{bn.fit}, \code{GBN} or \code{CI}.
#' @param ... parameters specific to the class used.
#'
#' @seealso \code{\link{KL.GBN}}, \code{\link{KL.CI}}, \code{\link{Fro.CI}}, \code{\link{Fro.GBN}}, \code{\link{Jeffreys.GBN}}, \code{\link{Jeffreys.CI}}
#'
#'@return A dataframe whose columns depend of the class of the object.
#'@export
#'

Jeffreys <- function (x, ...) {
  UseMethod("Jeffreys", x)
}


#' Jeffreys Divergence for \code{GBN}
#'
#'  \code{Jeffreys.GBN} returns the Jeffreys divergence between an object of class \code{GBN}  and its update after a standard parameter variation.
#'
#' Computation of the Jeffreys divergence between a Bayesian network and the additively perturbed Bayesian network, where the perturbation is either to the mean vector or to the covariance matrix.
#'
#'@return A dataframe including in the first column the variations performed and in the second column the corresponding Jeffreys divergences.
#'
#'@param x object of class \code{GBN}.
#'@param where character string: either \code{mean} or \code{covariance} for variations of the mean vector and covariance matrix respectively.
#'@param entry if \code{where == "mean"}, \code{entry} is the index of the entry of the mean vector to vary. If \code{where == "covariance"}, entry is a vector of length 2 indicating the entry of the covariance matrix to vary.
#'@param delta numeric vector, including the variation parameters that act additively.
#'@param ... additional arguments for compatibility.
#'
#'@references Goergen, C., & Leonelli, M. (2018). Model-preserving sensitivity analysis for families of Gaussian distributions. arXiv preprint arXiv:1809.10794.
#'
#'@seealso \code{\link{KL.GBN}}\code{\link{KL.CI}}, \code{\link{Fro.CI}}, \code{\link{Fro.GBN}}, \code{\link{Jeffreys.CI}}
#'
#'@examples Jeffreys(synthetic_gbn,"mean",2,seq(-1,1,0.1))
#'@examples Jeffreys(synthetic_gbn,"covariance",c(3,3),seq(-1,1,0.1))
#'
#'@importFrom ggplot2 ggplot
#'@importFrom ggplot2 geom_line
#'@importFrom ggplot2 geom_point
#'@importFrom ggplot2 labs
#'@importFrom ggplot2 aes
#'@export
#'
#'
Jeffreys.GBN <- function(x,where,entry,delta,  ...){
  gbn <- x
  if(where != "mean" & where!= "covariance") stop("where is either mean or covariance")
  KL <- numeric(length(delta))
  if(where == "mean"){
    d<-rep(0,length(gbn$mean))
    for(i in 1:length(KL)){
      d[entry] <- delta[i]
      KL[i]<- t(d)%*%solve(gbn$covariance)%*%d
    }
  }
  if(where == "covariance"){
    D <- matrix(0,length(gbn$mean),length(gbn$mean))
    inv_or <- solve(gbn$covariance)
    det_or <- det(gbn$covariance)
    for(i in 1:length(KL)){
      D[entry[1],entry[2]]<- delta[i]
      D[entry[2],entry[1]]<- delta[i]
      det_new <- det(gbn$covariance + D)
      if(is.psd(round(gbn$covariance+D,5)) & det_new > 1e-10){
        KL[i] <- 0.5*(sum(diag(inv_or%*%D))+log(det_or/det_new)) + 0.5*(sum(diag(solve(gbn$covariance + D)%*%gbn$covariance)) - nrow(gbn$covariance) + log(det(gbn$covariance + D)/det_or))
      }
      else{KL[i]<-NA}
    }
  }
  KL <- data.frame(Variation = delta,Jeffreys=KL)
  Variation <- KL$Variation
  Jeffreys <- KL$Jeffreys
    if(nrow(KL)==1){
      plot <- ggplot(data = KL, mapping = aes(x = Variation, y = Jeffreys)) + geom_point( na.rm = T) + labs(x = "delta",  y = "Jeffreys", title = "Jeffreys divergence") + theme_minimal()
    }else{
      plot <- ggplot(data = KL, mapping = aes(x = Variation, y = Jeffreys)) + geom_line( na.rm = T) + labs(x = "delta",  y = "Jeffreys", title = "Jeffreys divergence") + theme_minimal()
    }
  out <- list(Jeffreys = KL, plot = plot)
  attr(out,'class') <- 'jeffreys'
  return(out)
}

#' Jeffreys Divergence for \code{CI}
#'
#' \code{Jeffreys.CI} returns the Jeffreys divergence between an object of class \code{CI}  and its update after a model-preserving parameter variation.
#'
#' Computation of the Jeffreys divergence between a Bayesian network and its updated version after a model-preserving variation.
#'
#'@param x object of class \code{CI}.
#'@param type character string. Type of model-preserving co-variation: either \code{"total"}, \code{"partial"}, \code{row},\code{column} or \code{all}. If \code{all} the Jeffreys divergence is computed for every type of co-variation matrix.
#'@param entry a vector of length 2 indicating the entry of the covariance matrix to vary.
#'@param delta numeric vector with positive elements, including the variation parameters that act multiplicatively.
#'@param ... additional arguments for compatibility.
#'
#'@return A dataframe including in the first column the variations performed, and in the following columns the corresponding Jeffreys divergences for the chosen model-preserving co-variations.
#'
#'@references C. Görgen & M. Leonelli (2020), Model-preserving sensitivity analysis for families of Gaussian distributions.  Journal of Machine Learning Research, 21: 1-32.
#'
#'@seealso \code{\link{KL.GBN}}, \code{\link{KL.CI}}, \code{\link{Fro.CI}}, \code{\link{Fro.GBN}}, \code{\link{Jeffreys.GBN}}
#'
#'@examples Jeffreys(synthetic_ci,"total",c(1,1),seq(0.9,1.1,0.01))
#'@examples Jeffreys(synthetic_ci,"partial",c(1,4),seq(0.9,1.1,0.01))
#'@examples Jeffreys(synthetic_ci,"column",c(1,2),seq(0.9,1.1,0.01))
#'@examples Jeffreys(synthetic_ci,"row",c(3,2),seq(0.9,1.1,0.01))
#'@examples Jeffreys(synthetic_ci,"all",c(3,2),seq(0.9,1.1,0.01))
#'
#'@importFrom ggplot2 ggplot
#'@importFrom ggplot2 geom_line
#'@importFrom ggplot2 geom_point
#'@importFrom ggplot2 labs
#'@importFrom ggplot2 aes
#'@export
#'

Jeffreys.CI <- function(x, type, entry,delta, ...){
  ci <- x
  J <- numeric(length(delta))
  det_or <- det(ci$covariance)
  inv_or <- solve(ci$covariance)
  if(type == "total"){
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov <- total_covar_matrix(ci,entry,delta[i])
      det_new <- det(Cov*Delta*ci$covariance)
      if(is.psd(round(Cov*Delta*ci$covariance,2)) & det_new > 1e-10){
        J[i] <- 0.5*(log(det_or/det_new)-nrow(ci$covariance)+sum(diag(inv_or%*%(Cov*Delta*ci$covariance)))) + 0.5*(log(det_new/det_or)-nrow(ci$covariance)+ sum(diag(solve(Cov*Delta*ci$covariance)%*%ci$covariance)))
      }
      else{J[i]<- NA}
    }
  }
  if(type == "partial"){
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov <- partial_covar_matrix(ci,entry,delta[i])
      det_new <- det(Cov*Delta*ci$covariance)
      if(is.psd(round(Cov*Delta*ci$covariance,2))& det_new > 1e-10){
        J[i] <- 0.5*(log(det_or/det_new)-nrow(ci$covariance)+sum(diag(inv_or%*%(Cov*Delta*ci$covariance)))) + 0.5*(log(det_new/det_or)-nrow(ci$covariance)+ sum(diag(solve(Cov*Delta*ci$covariance)%*%ci$covariance)))
      }
      else{J[i]<- NA}
    }
  }
  if(type == "row"){
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov <- row_covar_matrix(ci, entry, delta[i])
      det_new <- det(Cov*Delta*ci$covariance)
      if(is.psd(round(Cov*Delta*ci$covariance,2))& det_new > 1e-10){
        J[i] <- 0.5*(log(det_or/det_new)-nrow(ci$covariance)+sum(diag(inv_or%*%(Cov*Delta*ci$covariance)))) + 0.5*(log(det_new/det_or)-nrow(ci$covariance)+ sum(diag(solve(Cov*Delta*ci$covariance)%*%ci$covariance)))
      }
      else{J[i]<- NA}
    }
  }
  if(type == "column"){
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov <- col_covar_matrix(ci,entry,delta[i])
      det_new <- det(Cov*Delta*ci$covariance)
      if(is.psd(round(Cov*Delta*ci$covariance,2))& det_new > 1e-10){
        J[i] <- 0.5*(log(det_or/det_new)-nrow(ci$covariance)+sum(diag(inv_or%*%(Cov*Delta*ci$covariance)))) + 0.5*(log(det_new/det_or)-nrow(ci$covariance)+ sum(diag(solve(Cov*Delta*ci$covariance)%*%ci$covariance)))
      }
      else{J[i]<- NA}
    }
  }
  if(type == "all"){
    J <- matrix(0,length(delta),4)
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov_col <- col_covar_matrix(ci,entry,delta[i])
      Cov_row <- row_covar_matrix(ci,entry,delta[i])
      Cov_par <- partial_covar_matrix(ci,entry,delta[i])
      Cov_tot <- total_covar_matrix(ci,entry,delta[i])
      det_new_col <- det(Cov_col*Delta*ci$covariance)
      det_new_row <- det(Cov_row*Delta*ci$covariance)
      det_new_par <- det(Cov_par*Delta*ci$covariance)
      det_new_tot <- det(Cov_tot*Delta*ci$covariance)
      if(is.psd(round(Cov_tot*Delta*ci$covariance,2))& det_new_tot > 1e-10){
        J[i,1] <- 0.5*(log(det_or/det_new_tot)-nrow(ci$covariance)+sum(diag(inv_or%*%(Cov_tot*Delta*ci$covariance)))) + 0.5*(log(det_new_tot/det_or)-nrow(ci$covariance)+ sum(diag(solve(Cov_tot*Delta*ci$covariance)%*%ci$covariance)))
      }
      else{J[i,1]<- NA}
      if(is.psd(round(Cov_par*Delta*ci$covariance,2))& det_new_par > 1e-10){
        J[i,2] <- 0.5*(log(det_or/det_new_par)-nrow(ci$covariance)+sum(diag(inv_or%*%(Cov_par*Delta*ci$covariance)))) + 0.5*(log(det_new_par/det_or)-nrow(ci$covariance)+ sum(diag(solve(Cov_par*Delta*ci$covariance)%*%ci$covariance)))
      }
      else{J[i,2]<- NA}
      if(is.psd(round(Cov_row*Delta*ci$covariance,2))& det_new_row > 1e-10){
        J[i,3] <- 0.5*(log(det_or/det_new_row)-nrow(ci$covariance)+sum(diag(inv_or%*%(Cov_row*Delta*ci$covariance)))) + 0.5*(log(det_new_row/det_or)-nrow(ci$covariance)+ sum(diag(solve(Cov_row*Delta*ci$covariance)%*%ci$covariance)))
      }
      else{J[i,3]<- NA}
      if(is.psd(round(Cov_col*Delta*ci$covariance,2))& det_new_col > 1e-10){
        J[i,4] <- 0.5*(log(det_or/det_new_col)-nrow(ci$covariance)+sum(diag(inv_or%*%(Cov_col*Delta*ci$covariance)))) + 0.5*(log(det_new_col/det_or)-nrow(ci$covariance)+ sum(diag(solve(Cov_col*Delta*ci$covariance)%*%ci$covariance)))
      }
      else{J[i,4]<- NA}
    }
  }
  if(type == "all"){
    J_data <- data.frame(Variation = delta, Total = J[,1], Partial = J[,2], Row_based = J[,3], Column_based = J[,4])}
   else{
    J_data <- data.frame(Variation = delta, Jeffreys= J)
    Variation <- J_data$Variation
    Jeffreys <- J_data$Jeffreys}
    if(type == "all"){
      ci <- gather(J_data, key = "scheme", value = "value", - Variation)
      Variation <- ci$Variation
      scheme <- ci$scheme
      value <- ci$value
      if(nrow(J_data) == 1){
        plot <- ggplot(data = ci, mapping = aes(x = Variation, y = value)) + geom_point(aes(color = scheme)) + labs( x = "delta", y = "Jeffreys", title = "Jeffreys divergence") + theme_minimal()
      } else{
        plot <- ggplot(data = ci, mapping = aes(x = Variation, y = value)) + geom_line(aes(color = scheme)) + labs( x = "delta", y = "Jeffreys", title = "Jeffreys divergence") + theme_minimal()
      }
    } else{
      if(nrow(J_data) == 1){
        plot <- ggplot(data = J_data, mapping = aes(x = Variation, y = Jeffreys)) + geom_point( na.rm = T) + labs(x = "delta",  y = "Jeffreys", title = "Jeffreys divergence") + theme_minimal()
      }else{
        plot <- ggplot(data = J_data, mapping = aes(x = Variation, y = Jeffreys)) + geom_line(na.rm = T ) + labs(x = "delta",  y = "Jeffreys", title = "Jeffreys divergence") + theme_minimal()
      }
    }
  out <- list(Jeffreys = J_data, plot = plot)
  attr(out,'class') <- 'jeffreys'
  return(out)
}


