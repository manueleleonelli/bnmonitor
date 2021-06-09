#' Bounds for the KL-divergence
#'
#' Computation of the bounds of the KL-divergence for variations of each parameter of a \code{CI} object.
#'
#'
#' Let \eqn{\Sigma} be the covariance matrix of a Gaussian Bayesian network with \eqn{n} vertices.
#' Let \eqn{D} and \eqn{\Delta} be  variation matrices acting additively on \eqn{\Sigma}. Let also \eqn{\tilde\Delta} be a model-preserving co-variation matrix.
#' Denote with \eqn{Y} and \eqn{\tilde{Y}} the original and the perturbed random vectors. Then for a standard sensitivity analysis
#' \deqn{KL(\tilde{Y}||Y)\leq 0.5n\max\left\{f(\lambda_{\max}(D\Sigma^{-1})),f(\lambda_{\min}(D\Sigma^{-1}))\right\}}
#' whilst for a model-preserving one
#' \deqn{KL(\tilde{Y}||Y)\leq 0.5n\max\left\{f(\lambda_{\max}(\tilde\Delta\circ\Delta)),f(\lambda_{\min}(\tilde\Delta\circ\Delta))\right\}}
#' where \eqn{\lambda_{\max}} and \eqn{\lambda_{\min}} are the largest and the smallest eigenvalues, respectively, \eqn{f(x)=\ln(1+x)-x/(1+x)} and \eqn{\circ} denotes the Schur or element-wise product.
#'@seealso \code{\link{KL.CI}}, \code{\link{KL.CI}}
#'
#'
#'@examples KL_bounds(synthetic_ci,1.05)
#'
#'
#'@param ci object of class \code{CI}.
#'@param delta multiplicative variation coefficient for the entry of the covariance matrix given in \code{entry}.
#'
#'@return A dataframe including the KL-divergence bound for each co-variation scheme (model-preserving and standard) and every entry of the covariance matrix. For variations leading to non-positive semidefinite matrix, the dataframe includes a \code{NA}.
#'
#'@references C. Görgen & M. Leonelli (2020), Model-preserving sensitivity analysis for families of Gaussian distributions.  Journal of Machine Learning Research, 21: 1-32.
#'@export

KL_bounds <- function(ci,delta){
  row <- c()
  column <- c()
  partial <- c()
  rowb <- c()
  colb <- c()
  standard <- c()
  total_mat <- matrix(delta,nrow=length(ci$mean),ncol=length(ci$mean))
  eig_total <- c(max(eigen(total_mat)$values),min(eigen(total_mat)$values))
  total <- 0.5*length(ci$mean)*rep(max(log(1+eig_total[1])-eig_total[1]/(1+eig_total[1]),log(1+eig_total[2])-eig_total[2]/(1+eig_total[2])),length(ci$mean)*(length(ci$mean)+1)/2)
  for(i in 1:length(ci$mean)){
    for(j in i:length(ci$mean)){
      d <- (delta -1)*ci$covariance[i,j]
      D <- matrix(0,nrow=length(ci$mean),ncol=length(ci$mean))
      D[i,j] <- d
      D[j,i] <- d
      row <- c(row,i)
      column <- c(column,j)
      partial_mat <- partial_covar_matrix(ci,c(i,j),delta)
      row_mat <- row_covar_matrix(ci,c(i,j),delta)
      column_mat <- row_covar_matrix(ci,c(i,j),delta)
      partial_mat[i,j] <- delta
      partial_mat[j,i] <- delta
      row_mat[i,j] <- delta
      row_mat[j,i] <- delta
      column_mat[i,j] <- delta
      column_mat[j,i] <- delta
      if(is.psd(round(D+ci$covariance,5))){
      eig_s <- c(max(eigen(D%*%solve(ci$covariance))$values),min(eigen(D%*%solve(ci$covariance))$values))
      standard <- c(standard,0.5*length(ci$mean)*max(log(1+eig_s[1])-eig_s[1]/(1+eig_s[1]),log(1+eig_s[2])-eig_s[2]/(1+eig_s[2])))
      }
      else{standard <- c(standard,NA)}
      if(is.psd(round(partial_mat*ci$covariance,5))){
      eig_partial <- c(max(eigen(partial_mat)$values),min(eigen(partial_mat)$values))
      partial <- c(partial, 0.5*length(ci$mean)*max(log(1+eig_partial[1])-eig_partial[1]/(1+eig_partial[1]),log(1+eig_partial[2])-eig_partial[2]/(1+eig_partial[2])))}
      else{partial <- c(partial,NA)}
      if(is.psd(round(row_mat*ci$covariance,5))){
      eig_row <- c(max(eigen(row_mat)$values),min(eigen(row_mat)$values))
      rowb <- c(rowb,0.5*length(ci$mean)*max(log(1+eig_row[1])-eig_row[1]/(1+eig_row[1]),log(1+eig_row[2])-eig_row[2]/(1+eig_row[2])))}
      else{rowb <- c(rowb,NA)}
      if(is.psd(round(row_mat*ci$covariance,5))){
      eig_col <- c(max(eigen(column_mat)$values),min(eigen(column_mat)$values))
      colb <- c(colb,0.5*length(ci$mean)*max(log(1+eig_col[1])-eig_col[1]/(1+eig_col[1]),log(1+eig_col[2])-eig_col[2]/(1+eig_col[2])))}
      else{colb <- c(colb,NA)}
    }
  }
  return(as.data.frame(cbind(row = row, col = column, standard = standard, total=total,partial=partial,row_based=rowb,col_based=colb)))
}

