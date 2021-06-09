#' KL Divergence
#'
#' \code{KL} returns the Kullback-Leibler (KL) divergence between a Bayesian network and its update after parameter variation.
#'
#' The details depend on the class the method \code{KL} is applied to.
#'
#' @param x object of class \code{bn.fit}, \code{GBN} or \code{CI}.
#' @param ... parameters specific to the class used.
#'
#' @seealso \code{\link{KL.GBN}}, \code{\link{KL.CI}}, \code{\link{Fro.CI}}, \code{\link{Fro.GBN}}, \code{\link{Jeffreys.GBN}}, \code{\link{Jeffreys.CI}}
#'
#'@return A dataframe whose columns depend of the class of the object.
#'@export
#'

KL <- function (x, ...) {
  UseMethod("KL", x)
}

#' KL Divergence for \code{bn.fit}
#'
#' \code{KL.bn.fit} returns the Kullback-Leibler (KL) divergence between a Bayesian network and its update after parameter variation.
#'
#'
#' The Bayesian network on which parameter variation is being conducted should be expressed as a \code{bn.fit} object.
#' The name of the node to be varied, its level and its parent's levels should be specified.
#' The parameter variation specified by the function is:
#'
#'  P ( \code{node} = \code{value_node} | parents = \code{value_parents} ) = \code{new_value}
#'
#'@seealso \code{\link{CD}}
#'
#'@param x object of class \code{bn.fit}.
#'@param node character string. Node of which the conditional probability distribution is being changed.
#'@param value_node character string. Level of \code{node}.
#'@param value_parents character string. Levels of \code{node}'s parents. The levels should be defined according to the order of the parents in \code{bnfit[[node]][["parents"]]}. If \code{node} has no parents, then it should be set to \code{NULL}.
#'@param new_value numeric vector with elements between 0 and 1. Values to which the parameter should be updated. It can take a specific value or more than one. In the case of more than one value, these should be defined through a vector with an increasing order of the elements. \code{new_value} can also be set to the character string \code{all}: in this case a sequence of possible parameter changes ranging from 0.05 to 0.95 is considered.
#'@param covariation character string. Co-variation scheme to be used for the updated Bayesian network. Can take values \code{uniform}, \code{proportional}, \code{orderp}, \code{all}. If equal to \code{all}, uniform, proportional and order-preserving co-variation schemes are used. Set by default to \code{proportional}.
#'@param ... additional parameters to be added to the plot.
#'
#'
#' @return A dataframe with the varied parameter and the KL divergence for different co-variation schemes. If \code{plot = TRUE} the function returns a plot of the KL divergences.
#'
#'@references Kullback, S., & Leibler, R. A. (1951). On information and sufficiency. The annals of mathematical statistics, 22(1), 79-86.
#'@references Leonelli, M., Goergen, C., & Smith, J. Q. (2017). Sensitivity analysis in multilinear probabilistic models. Information Sciences, 411, 84-97.
#'
#'@examples KL(synthetic_bn, "y2", "1", "2", "all", "all")
#'@examples KL(synthetic_bn, "y1", "2", NULL, 0.3, "all")
#'
#'@importClassesFrom bnlearn bn.fit
#'@importFrom stats coef
#'@importFrom bnlearn as.grain
#'@importFrom gRain querygrain
#'@importFrom graphics lines points
#'@importFrom gRbase compile
#'@importFrom ggplot2 ggplot
#'@importFrom ggplot2 geom_line
#'@importFrom ggplot2 geom_point
#'@importFrom ggplot2 labs
#'@importFrom ggplot2 aes
#'@export
#'

KL.bn.fit <-
  function(x,
           node,
           value_node,
           value_parents,
           new_value,
           covariation = "proportional",
           ...) {
    bnfit <- x
    suppressWarnings(if (new_value == "all") {
      new_value2 <-
        sort(c(seq(0.05, 0.95, by = 0.05), coef(bnfit[[node]])[t(append(value_node, value_parents))]))
    } else{
      new_value2 <- new_value
    })
    new_value_op <- rep(NA, length(new_value2))
    bnfit.new <- vector("list", length = length(new_value2))
    if (covariation == "all") {
      bnfit.new.scheme <- vector("list", length = 3)
    } else{
      bnfit.new.scheme <- vector("list", length = 1)
    }
    if (covariation == "uniform" || covariation == "all") {
      for (i in 1:length(new_value2)) {
        bnfit.new[[i]] <- uniform_covar(
          bnfit = bnfit,
          node = node,
          value_node = value_node,
          value_parents = value_parents,
          new_value = new_value2[i]
        )
      }
      bnfit.new.scheme[[1]] <- bnfit.new

    }
    if (covariation == "proportional" || covariation == "all") {
      for (i in 1:length(new_value2)) {
        bnfit.new[[i]] <- proportional_covar(
          bnfit = bnfit,
          node = node,
          value_node = value_node,
          value_parents = value_parents,
          new_value = new_value2[i]
        )
      }
      if (covariation == "all") {
        bnfit.new.scheme[[2]] <- bnfit.new
      } else{
        bnfit.new.scheme[[1]] <- bnfit.new
      }
    }
    if (covariation == "orderp" || covariation == "all") {
      new_value_op <- new_value2
      length.node <- dim(coef(bnfit[[node]]))[1]
      node_values <- numeric(length.node)
      cpt.new <- coef(bnfit[[node]])
      for (i in 1:length.node) {
        node_values[i] <-
          cpt.new[t(append(unlist(dimnames(coef(
            bnfit[[node]]
          ))[1])[i], value_parents))]
      }
      order <- rank(node_values, ties.method = "first")
      s <-
        order[unname(which(unlist(dimnames(coef(
          bnfit[[node]]
        ))[1]) == value_node))]
      if (s == length.node & covariation == "orderp") {
        stop (
          "Order preserving covariation scheme can't work by varying the last parameter in the increasing order of the parameters."
        )
      }
      if (s == length.node & covariation == "all") {
        warning(
          "The KL-divergence with order-preserving scheme won't be returned since the last parameter in the increasing order of the parameters is being varied."
        )
        bnfit.new.scheme[[3]] <- NULL
      }
      if (s != length.node) {
        new_value_op <-
          new_value_op[new_value_op <= 1 / (1 + dim(coef(bnfit[[node]]))[1] - s)]
        if (length(new_value_op) != 0) {
          bnfit.new.op <- vector("list", length = length(new_value_op))
          for (i in 1:length(new_value_op)) {
            bnfit.new.op[[i]] <- orderp_covar(
              bnfit = bnfit,
              node = node,
              value_node = value_node,
              value_parents = value_parents,
              new_value = new_value_op[i]
            )
          }
          if (covariation == "all") {
            bnfit.new.scheme[[3]] <- bnfit.new.op
          } else{
            bnfit.new.scheme[[1]] <- bnfit.new.op
          }
        }
      }
    }
    if (covariation == "all") {
      if (length(bnfit.new.scheme) == 2) {
        KL <-
          data.frame(new_value2, rep(0, length(new_value2)), rep(0, length(new_value2)))
        colnames(KL) <- c("New_value", "Uniform", "Proportional")
      } else{
        KL <-
          data.frame(new_value2,
                     rep(0, length(new_value2)),
                     rep(0, length(new_value2)),
                     c(rep(0, length(new_value_op)), rep(
                       NA, length(new_value2) - length(new_value_op)
                     )))
        colnames(KL) <-
          c("New_value",
            "Uniform",
            "Proportional",
            "Order Preserving")
      }
    } else{
      if (covariation == "orderp") {
        KL <-
          data.frame(new_value2, c(rep(0, length(new_value_op)), rep(
            NA, length(new_value2) - length(new_value_op)
          )))
        colnames(KL) <- c("New_value", "KL_DIVERGENCE")
      } else{
        KL <- data.frame(new_value2, rep(0, length(new_value2)))
        colnames(KL) <- c("New_value", "KL_DIVERGENCE")
      }
    }
    grain.bn <- compile(as.grain(bnfit))
    if (is.null(value_parents)) {
      pr_parents <- 1
    } else{
      pr_parents <-
        querygrain(grain.bn, nodes = as.vector(bnfit[[node]][["parents"]]), type =
                     "joint")[t(value_parents)]
    }
    if (length(new_value2) == length(new_value_op)) {
      for (k in 1:length(bnfit.new.scheme)) {
        for (j in 1:length(new_value2)) {
          for (i in 1:dim(coef(bnfit[[node]]))[1]) {
            KL[j, k + 1] <-
              KL[j, k + 1] + (coef(bnfit[[node]])[t(append(unlist(dimnames(
                coef(bnfit[[node]])
              )[1])[i], value_parents))]) * log((coef(bnfit[[node]])[t(append(unlist(dimnames(
                coef(bnfit[[node]])
              )[1])[i], value_parents))]) / coef(bnfit.new.scheme[[k]][[j]][[node]])[t(append(unlist(dimnames(
                coef(bnfit.new.scheme[[k]][[j]][[node]])
              )[1])[i], value_parents))])
          }
          KL[j, k + 1] <- pr_parents * KL[j, k + 1]
        }
      }
    } else{
      if (covariation == "all") {
        for (k in 1:2) {
          for (j in 1:length(new_value2)) {
            quotients <- c()
            for (i in 1:dim(coef(bnfit[[node]]))[1]) {
              KL[j, k + 1] <-
                KL[j, k + 1] + (coef(bnfit[[node]])[t(append(unlist(dimnames(
                  coef(bnfit[[node]])
                )[1])[i], value_parents))]) * log((coef(bnfit[[node]])[t(append(unlist(dimnames(
                  coef(bnfit[[node]])
                )[1])[i], value_parents))]) / coef(bnfit.new.scheme[[k]][[j]][[node]])[t(append(unlist(dimnames(
                  coef(bnfit.new.scheme[[k]][[j]][[node]])
                )[1])[i], value_parents))])
            }
            KL[j, k + 1] <- pr_parents * KL[j, k + 1]
          }
        }
        if (length(new_value_op) != 0) {
          for (t in 1:(length(new_value_op))) {
            quotients <- c()
            for (i in 1:dim(coef(bnfit[[node]]))[1]) {
              KL[t, 4] <-
                KL[t, 4] + (coef(bnfit[[node]])[t(append(unlist(dimnames(
                  coef(bnfit[[node]])
                )[1])[i], value_parents))]) * log((coef(bnfit[[node]])[t(append(unlist(dimnames(
                  coef(bnfit[[node]])
                )[1])[i], value_parents))]) / coef(bnfit.new.scheme[[3]][[t]][[node]])[t(append(unlist(dimnames(
                  coef(bnfit.new.scheme[[3]][[t]][[node]])
                )[1])[i], value_parents))])
            }
            KL[t, 4] <- pr_parents * KL[t, 4]
          }
        }
      } else{
        if (length(new_value_op) != 0) {
          for (t in 1:(length(new_value_op))) {
            quotients <- c()
            for (i in 1:dim(coef(bnfit[[node]]))[1]) {
              KL[t, 2] <-
                KL[t, 2] + (coef(bnfit[[node]])[t(append(unlist(dimnames(
                  coef(bnfit[[node]])
                )[1])[i], value_parents))]) * log((coef(bnfit[[node]])[t(append(unlist(dimnames(
                  coef(bnfit[[node]])
                )[1])[i], value_parents))]) / coef(bnfit.new.scheme[[1]][[t]][[node]])[t(append(unlist(dimnames(
                  coef(bnfit.new.scheme[[1]][[t]][[node]])
                )[1])[i], value_parents))])
            }
            KL[t, 2] <- pr_parents * KL[t, 2]
          }
        }
      }
    }
   # if (all(is.na(KL[,-1]))) {
  #    plot <- FALSE
  #    warning("The plot won't be showed since all the values are not possible")
   # }
    ci <- gather(KL, key = "scheme", value = "value", - New_value)
    New_value <- ci$New_value
    scheme <- ci$scheme
    value <- ci$value
    if (covariation == "all") {
      if (nrow(KL) == 1) {
        if(ncol(KL) == 3){
          plot <- ggplot(data = ci, mapping = aes(x = New_value, y = value)) + geom_point(aes(color = scheme))  + labs(x = "new value", y = "CD", title = "CD distance") + theme_minimal()
        } else{
          plot <- ggplot(data = ci, mapping = aes(x = New_value, y = value)) + geom_point(aes(color = scheme))  + labs(x = "new value", y = "CD", title = "CD distance") + theme_minimal()
        }
      }
      else{
        if(ncol(KL) == 3){
          plot <- ggplot(data = ci, mapping = aes(x = New_value, y = value)) + geom_line(aes(color = scheme))  + labs(x = "new value", y = "CD", title = "CD distance") + theme_minimal()
        } else{
          plot <- ggplot(data = ci, mapping = aes(x = New_value, y = value)) + geom_line(aes(color = scheme))  + labs(x = "new value", y = "CD", title = "CD distance") + theme_minimal()
        }
      }
    }
       else{
        if (nrow(KL) == 1) {
          plot <- ggplot(data = KL, mapping = aes(x = KL[,1], y = KL[,2])) + geom_point( na.rm = T) + labs(x = "New value", y = "KL", title = "KL divergence") + theme_minimal()
        } else{
          plot <- ggplot(data = KL, mapping = aes(x = KL[,1], y = KL[,2])) + geom_line( na.rm = T) + labs(x = "New value", y = "KL", title = "KL divergence") + theme_minimal()
        }
       }
    out <- list(KL = KL, plot = plot)
    attr(out,'class') <- 'kl'
    return(out)
  }






#' KL Divergence for \code{CI}
#'
#' \code{KL.CI} returns the Kullback-Leibler (KL) divergence between an object of class \code{CI}  and its update after a model-preserving parameter variation.
#'
#' Computation of the KL divergence between a Bayesian network and its updated version after a model-preserving variation.
#'
#'@param x object of class \code{CI}.
#'@param type character string. Type of model-preserving co-variation: either \code{"total"}, \code{"partial"}, \code{row},\code{column} or \code{all}. If \code{all} the KL divergence is computed for every type of co-variation matrix.
#'@param entry a vector of length 2 indicating the entry of the covariance matrix to vary.
#'@param delta numeric vector with positive elements, including the variation parameters that act multiplicatively.
#'@param ... additional arguments for compatibility.
#'
#'@return A dataframe including in the first column the variations performed, and in the following columns the corresponding KL divergences for the chosen model-preserving co-variations.
#'
#'@references C. Görgen & M. Leonelli (2020), Model-preserving sensitivity analysis for families of Gaussian distributions.  Journal of Machine Learning Research, 21: 1-32.
#'@seealso \code{\link{KL.GBN}}, \code{\link{Fro.CI}}, \code{\link{Fro.GBN}}, \code{\link{Jeffreys.GBN}}, \code{\link{Jeffreys.CI}}
#'
#'@examples KL(synthetic_ci, "total", c(1,1), seq(0.9,1.1,0.01))
#'@examples KL(synthetic_ci, "partial", c(1,4), seq(0.9,1.1,0.01))
#'@examples KL(synthetic_ci, "column", c(1,2), seq(0.9,1.1,0.01))
#'@examples KL(synthetic_ci, "row", c(3,2), seq(0.9,1.1,0.01))
#'@examples KL(synthetic_ci, "all", c(3,2), seq(0.9,1.1,0.01))
#'
#'@importFrom ggplot2 ggplot
#'@importFrom ggplot2 geom_line
#'@importFrom ggplot2 geom_point
#'@importFrom ggplot2 labs
#'@importFrom ggplot2 aes
#'@export
#'

KL.CI <- function(x, type, entry, delta,  ...){
  ci <- x
  KL <- numeric(length(delta))
  det_or <- det(ci$covariance)
  inv_or <- solve(ci$covariance)
  if(type == "total"){
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov <- total_covar_matrix(ci,entry,delta[i])
      det_new <- det(Cov*Delta*ci$covariance)
      if(is.psd(round(Cov*Delta*ci$covariance,2))& det_new > 1e-10){
        KL[i] <- 0.5*(log(det_or/det(Cov*Delta*ci$covariance))-nrow(ci$covariance)+sum(diag(inv_or%*%(Cov*Delta*ci$covariance))))
      }
      else{KL[i]<- NA}
    }
  }
  if(type == "partial"){
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov <- partial_covar_matrix(ci,entry,delta[i])
      det_new <- det(Cov*Delta*ci$covariance)
      if(is.psd(round(Cov*Delta*ci$covariance,2))& det_new > 1e-10){
        KL[i] <- 0.5*(log(det_or/det(Cov*Delta*ci$covariance))-nrow(ci$covariance)+sum(diag(inv_or%*%(Cov*Delta*ci$covariance))))
      }
      else{KL[i]<- NA}
    }
  }
  if(type == "row"){
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov <- row_covar_matrix(ci, entry, delta[i])
      det_new <- det(Cov*Delta*ci$covariance)
      if(is.psd(round(Cov*Delta*ci$covariance,2))& det_new > 1e-10){
        KL[i] <- 0.5*(log(det_or/det(Cov*Delta*ci$covariance))-nrow(ci$covariance)+sum(diag(inv_or%*%(Cov*Delta*ci$covariance))))
      }
      else{KL[i]<- NA}
    }
   }
  if(type == "column"){
    for(i in 1:length(delta)){
      Delta <- variation_mat(ci,entry,delta[i])
      Cov <- col_covar_matrix(ci,entry,delta[i])
      det_new <- det(Cov*Delta*ci$covariance)
      if(is.psd(round(Cov*Delta*ci$covariance,2))& det_new > 1e-10){
        KL[i] <- 0.5*(log(det_or/det(Cov*Delta*ci$covariance))-nrow(ci$covariance)+sum(diag(inv_or%*%(Cov*Delta*ci$covariance))))
      }
      else{KL[i]<- NA}
    }
  }
  if(type == "all"){
    KL <- matrix(0,length(delta),4)
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
      if(is.psd(round(Cov_tot*Delta*ci$covariance,2))&det_new_tot > 1e-10){
        KL[i,1] <- 0.5*(log(det_or/det(Cov_tot*Delta*ci$covariance))-nrow(ci$covariance)+sum(diag(inv_or%*%(Cov_tot*Delta*ci$covariance))))
      }
      else{KL[i,1]<- NA}
      if(is.psd(round(Cov_par*Delta*ci$covariance,2))&det_new_par > 1e-10){
        KL[i,2] <- 0.5*(log(det_or/det(Cov_par*Delta*ci$covariance))-nrow(ci$covariance)+sum(diag(inv_or%*%(Cov_par*Delta*ci$covariance))))
      }
      else{KL[i,2]<- NA}
      if(is.psd(round(Cov_row*Delta*ci$covariance,2))&det_new_row > 1e-10){
        KL[i,3] <- 0.5*(log(det_or/det(Cov_row*Delta*ci$covariance))-nrow(ci$covariance)+sum(diag(inv_or%*%(Cov_row*Delta*ci$covariance))))
      }
      else{KL[i,3]<- NA}
      if(is.psd(round(Cov_col*Delta*ci$covariance,2))&det_new_col > 1e-10){
        KL[i,4] <- 0.5*(log(det_or/det(Cov_col*Delta*ci$covariance))-nrow(ci$covariance)+sum(diag(inv_or%*%(Cov_col*Delta*ci$covariance))))
      }
      else{KL[i,4]<- NA}
    }
  }
  if(type == "all"){KL_data <- data.frame(Variation = delta, Total = KL[,1], Partial = KL[,2], Row_based = KL[,3], Column_based = KL[,4])}
  else{ KL_data <- data.frame(Variation = delta, KL= KL)
  KL <- KL_data$KL}
  ci <- gather(KL_data, key = "scheme", value = "value", - Variation)
  Variation <- ci$Variation
  scheme <- ci$scheme
  value <- ci$value
    if(type == "all"){
      if(nrow(KL_data) == 1){
        plot <- ggplot(data = ci, mapping = aes(x = Variation, y = value)) + geom_point(aes(color = scheme))  + labs(x = "delta", y = "KL", title = "KL divergence") + theme_minimal()
      } else{
        plot <- ggplot(data = ci, mapping = aes(x = Variation, y = value)) + geom_line(aes(color = scheme)) + labs(x = "delta", y = "KL", title = "KL divergence") + theme_minimal()
      }
    } else{
      if(nrow(KL_data) == 1){
        plot <- ggplot(data = KL_data, mapping = aes(x = Variation, y = KL)) + geom_point( na.rm = T) + labs(x = "delta", y = "KL", title = "KL divergence") + theme_minimal()
      }else{
       plot <- ggplot(data = KL_data, mapping = aes(x = Variation, y = KL)) + geom_line(na.rm = T ) + labs(x = "delta", y = "KL", title = "KL divergence") + theme_minimal()
      }
    }
  out <- list(KL = KL_data, plot = plot)
  attr(out,'class') <- 'kl'
  return(out)
}


#' KL Divergence for \code{GBN}
#'
#'  \code{KL.GBN} returns the Kullback-Leibler (KL) divergence between an object of class \code{GBN}  and its update after a standard parameter variation.
#'
#' Computation of the KL divergence between a Bayesian network and the additively perturbed Bayesian network, where the perturbation is either to the mean vector or to the covariance matrix.
#'
#'@return A dataframe including in the first column the variations performed and in the second column the corresponding KL divergences.
#'
#'@param x object of class \code{GBN}.
#'@param where character string: either \code{mean} or \code{covariance} for variations of the mean vector and covariance matrix respectively.
#'@param entry if \code{where == "mean"}, \code{entry} is the index of the entry of the mean vector to vary. If \code{where == "covariance"}, entry is a vector of length 2 indicating the entry of the covariance matrix to vary.
#'@param delta numeric vector, including the variation parameters that act additively.
#'@param ... additional arguments for compatibility.
#'
#'@references Gómez-Villegas, M. A., Maín, P., & Susi, R. (2007). Sensitivity analysis in Gaussian Bayesian networks using a divergence measure. Communications in Statistics—Theory and Methods, 36(3), 523-539.
#'@references Gómez-Villegas, M. A., Main, P., & Susi, R. (2013). The effect of block parameter perturbations in Gaussian Bayesian networks: Sensitivity and robustness. Information Sciences, 222, 439-458.
#'
#'@seealso \code{\link{KL.CI}}, \code{\link{Fro.CI}}, \code{\link{Fro.GBN}}, \code{\link{Jeffreys.GBN}}, \code{\link{Jeffreys.CI}}
#'
#'@examples KL(synthetic_gbn,"mean",2,seq(-1,1,0.1))
#'@examples KL(synthetic_gbn,"covariance",c(3,3),seq(-1,1,0.1))
#'
#'@importFrom ggplot2 ggplot
#'@importFrom ggplot2 geom_line
#'@importFrom ggplot2 geom_point
#'@importFrom ggplot2 labs
#'@importFrom ggplot2 aes
#'@export
#'
#'
KL.GBN <- function(x,where,entry,delta,  ...){
  gbn <- x
  if(where != "mean" & where!= "covariance") stop("where is either mean or covariance")
  KL <- numeric(length(delta))
  if(where == "mean"){
    inv_or <- solve(gbn$covariance)
    d<-rep(0,length(gbn$mean))
    for(i in 1:length(KL)){
      d[entry] <- delta[i]
      KL[i]<- 0.5*t(d)%*%inv_or%*%d
    }
  }
  if(where == "covariance"){
    D <- matrix(0,length(gbn$mean),length(gbn$mean))
    inv_or <- solve(gbn$covariance)
    det_or <- det(gbn$covariance)
    for(i in 1:length(KL)){
      D[entry[1],entry[2]]<- delta[i]
      D[entry[2],entry[1]]<- delta[i]
      if(is.psd(round(gbn$covariance+D),2) & det(gbn$covariance+D) > 1e-10){
        KL[i] <- 0.5*(sum(diag(inv_or%*%D))+log(det_or/det(gbn$covariance+D)))
      }
      else{KL[i]<-NA}
    }
  }
  KL_data <- data.frame(Variation = delta, KL=KL)
  Variation <- KL_data$Variation
  KL <- KL_data$KL
    if(nrow(KL_data)==1){
      plot <- ggplot(data = KL_data, mapping = aes(x = Variation, y = KL)) + geom_point( na.rm = T) + labs(x = "delta", y = "KL", title = "KL divergence") + theme_minimal()
    }else{
      plot <- ggplot(data = KL_data, mapping = aes(x = Variation, y = KL)) + geom_line( na.rm = T) + labs(x = "delta", y = "KL", title = "KL divergence") + theme_minimal()
    }
  out <- list(KL = KL_data, plot = plot)
  attr(out,'class') <- 'kl'
  return(out)
}


