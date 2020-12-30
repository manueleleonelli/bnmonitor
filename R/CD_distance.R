#' CD-distance
#'
#' Chan-Darwiche (CD) distance between a Bayesian network and its update after parameter variation.
#'
#' The Bayesian network on which parameter variation is being conducted should be expressed as a \code{bn.fit} object.
#' The name of the node to be varied, its level and its parent's levels should be specified.
#' The parameter variation specified by the function is:
#'
#'  P ( \code{node} = \code{value_node} | parents = \code{value_parents} ) = \code{new_value}
#'
#' The CD distance between two probability distributions \eqn{P} and \eqn{P'} defined over the same sample space \eqn{\mathcal{Y}} is defined as
#' \deqn{CD(P,P')= \log\max_{y\in\mathcal{Y}}\left(\frac{P(y)}{P'(y)}\right) - \log\min_{y\in\mathcal{Y}}\left(\frac{P(y)}{P'(y)}\right)}
#'
#' @seealso \code{\link{KL.bn.fit}}
#'
#'@param bnfit object of class \code{bn.fit}.
#'@param node character string. Node of which the conditional probability distribution is being changed.
#'@param value_node character string. Level of \code{node}.
#'@param value_parents character string. Levels of \code{node}'s parents. The levels should be defined according to the order of the parents in \code{bnfit[[node]][["parents"]]}. If \code{node} has no parents, then it should be set to \code{NULL}.
#'@param new_value numeric vector with elements between 0 and 1. Values to which the parameter should be updated. It can take a specific value or more than one. In the case of more than one value, these should be defined through a vector with an increasing order of the elements. \code{new_value} can also be set to the character string \code{all}: in this case a sequence of possible parameter changes ranging from 0.05 to 0.95 is considered.
#'@param covariation character string. Co-variation scheme to be used for the updated Bayesian network. Can take values \code{uniform}, \code{proportional}, \code{orderp}, \code{all}. If equal to \code{all}, uniform, proportional and order-preserving co-variation schemes are used. Set by default to \code{proportional}.
#'
#'@references Chan, H., & Darwiche, A. (2005). A distance measure for bounding probabilistic belief change. International Journal of Approximate Reasoning, 38(2), 149-174.
#'@references Renooij, S. (2014). Co-variation for sensitivity analysis in Bayesian networks: Properties, consequences and alternatives. International Journal of Approximate Reasoning, 55(4), 1022-1042.
#'
#'@return The function \code{CD} returns a dataframe including in the first column the variations performed, and in the following columns the corresponding CD distances for the chosen co-variation schemes.
#'
#'@examples CD(synthetic_bn, "y2", "1", "2", "all", "all")
#'@examples CD(synthetic_bn, "y1", "2", NULL, 0.3, "all")
#'
#'@importFrom stats coef
#'@importFrom graphics lines points
#'@importFrom tidyr gather
#'@importClassesFrom bnlearn bn.fit
#'@export
CD <- function(bnfit, node, value_node, value_parents, new_value, covariation = "proportional") {
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
          "The CD distance with order-preserving scheme won't be returned since the last parameter in the increasing order of the parameters is being varied."
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
        CD <-
          data.frame(new_value2, rep(NA, length(new_value2)), rep(NA, length(new_value2)))
        colnames(CD) <- c("New_value", "Uniform", "Proportional")
      } else{
        CD <-
          data.frame(new_value2,
                     rep(NA, length(new_value2)),
                     rep(NA, length(new_value2)),
                     rep(NA, length(new_value2)))
        colnames(CD) <-
          c("New_value",
            "Uniform",
            "Proportional",
            "Order Preserving")
      }
    } else{
      CD <- data.frame(new_value2, rep(NA, length(new_value2)))
      colnames(CD) <- c("New_value", "CD_distance")
    }

    if (length(new_value2) == length(new_value_op)) {
      for (k in 1:length(bnfit.new.scheme)) {
        for (j in 1:length(new_value2)) {
          quotients <- c()
          for (i in 1:dim(coef(bnfit[[node]]))[1]) {
            quotients[i] <-
              (coef(bnfit.new.scheme[[k]][[j]][[node]])[t(append(unlist(dimnames(
                coef(bnfit.new.scheme[[k]][[j]][[node]])
              )[1])[i], value_parents))]) / (coef(bnfit[[node]])[t(append(unlist(dimnames(
                coef(bnfit[[node]])
              )[1])[i], value_parents))])
          }
          CD[j, k + 1] <- log(max(quotients)) - log(min(quotients))
        }
      }
    } else{
      if (covariation == "all") {
        for (k in 1:2) {
          for (j in 1:length(new_value2)) {
            quotients <- c()
            for (i in 1:dim(coef(bnfit[[node]]))[1]) {
              quotients[i] <-
                (coef(bnfit.new.scheme[[k]][[j]][[node]])[t(append(unlist(dimnames(
                  coef(bnfit.new.scheme[[k]][[j]][[node]])
                )[1])[i], value_parents))]) / (coef(bnfit[[node]])[t(append(unlist(dimnames(
                  coef(bnfit[[node]])
                )[1])[i], value_parents))])
            }
            CD[j, k + 1] <-
              log(max(quotients)) - log(min(quotients))
          }
        }
        if (length(new_value_op) != 0) {
          for (t in 1:(length(new_value_op))) {
            quotients <- c()
            for (i in 1:dim(coef(bnfit[[node]]))[1]) {
              quotients[i] <-
                (coef(bnfit.new.scheme[[3]][[t]][[node]])[t(append(unlist(dimnames(
                  coef(bnfit.new.scheme[[3]][[t]][[node]])
                )[1])[i], value_parents))]) / (coef(bnfit[[node]])[t(append(unlist(dimnames(
                  coef(bnfit[[node]])
                )[1])[i], value_parents))])
            }
            CD[t, 4] <- log(max(quotients)) - log(min(quotients))
          }
        }
      } else{
        if (length(new_value_op) != 0) {
          for (t in 1:(length(new_value_op))) {
            quotients <- c()
            for (i in 1:dim(coef(bnfit[[node]]))[1]) {
              quotients[i] <-
                (coef(bnfit.new.scheme[[1]][[t]][[node]])[t(append(unlist(dimnames(
                  coef(bnfit.new.scheme[[1]][[t]][[node]])
                )[1])[i], value_parents))]) / (coef(bnfit[[node]])[t(append(unlist(dimnames(
                  coef(bnfit[[node]])
                )[1])[i], value_parents))])
            }
            CD[t, 2] <- log(max(quotients)) - log(min(quotients))
          }
        }
      }
    }
   # if (all(is.na(CD[,-1]))) {
  #    plot <- FALSE
  #    warning("The plot won't be showed since all the values are not possible")
  #  }
    ci <- gather(CD, key = "scheme", value = "value", - New_value)
    New_value <- ci$New_value
    scheme <- ci$scheme
    value <- ci$value
      if (covariation == "all") {
        if (nrow(CD) == 1) {
          if(ncol(CD) == 3){
            plot <- ggplot(data = ci, mapping = aes(x = New_value, y = value)) + geom_point(aes(color = scheme))  + labs(x = "new value", y = "CD", title = "CD distance") + theme_minimal()
          } else{
            plot <- ggplot(data = ci, mapping = aes(x = New_value, y = value)) + geom_point(aes(color = scheme))  + labs(x = "new value", y = "CD", title = "CD distance") + theme_minimal()
          }
        }
        else{
          if(ncol(CD) == 3){
            plot <- ggplot(data = ci, mapping = aes(x = New_value, y = value)) + geom_line(aes(color = scheme))  + labs(x = "new value", y = "CD", title = "CD distance") + theme_minimal()
          } else{
            plot <- ggplot(data = ci, mapping = aes(x = New_value, y = value)) + geom_line(aes(color = scheme))  + labs(x = "new value", y = "CD", title = "CD distance") + theme_minimal()
          }
        }
      }
      else{
        if (nrow(CD) == 1) {
          plot <- ggplot(data = CD, mapping = aes(x = CD[,1], y = CD[,2])) + geom_point( na.rm = T) + labs(x = "new value", y = "CD", title = "CD distance") + theme_minimal()
        } else{
          plot <- ggplot(data = CD, mapping = aes(x = CD[,1], y = CD[,2])) + geom_line( na.rm = T) + labs(x = "new value", y = "CD", title = "CD distance") + theme_minimal()
        }
      }
    out <- list(CD = CD, plot = plot)
    attr(out,'class') <- 'CD'
    return(out)
}
