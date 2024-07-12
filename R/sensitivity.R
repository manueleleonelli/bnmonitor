#' Sensitivity function
#'
#' \code{sensitivity} returns the sensitivity function for a probabilistic query of interest with respect to a parameter change defined by the user.
#'
#' The Bayesian network on which parameter variation is being conducted should be expressed as a bn.fit object.
#' The name of the node to be varied, its level and its parent's level should be specified.
#' The parameter variation specified by the function is:
#'
#'  P ( \code{node} = \code{value_node} | parents = \code{value_parents} ) = \code{new_value}
#'
#'and the probabilistic query of interest is:
#'
#'  P ( \code{interest_node} = \code{interest_node_value} | \code{evidence_nodes} = \code{evidence_states} )
#'
#' @seealso \code{\link{covariation}}, \code{\link{sensquery}}
#'
#' @return A dataframe with the varied parameter values and the output probabilities for the co-variation schemes selected. If \code{plot = TRUE} the function also returns a plot of the sensitivity function.
#'
#'@param bnfit object of class \code{bn.fit}.
#'@param interest_node character string. Node of the probability query of interest.
#'@param interest_node_value character string. Level of \code{interest_node}.
#'@param evidence_nodes character string. Evidence nodes. If \code{NULL} no evidence is considered. Set by default to \code{NULL}.
#'@param evidence_states character string. Levels of \code{evidence_nodes}. If \code{NULL} no evidence is considered. If \code{evidence_nodes="NULL"}, \code{evidence_states} should be set to \code{NULL}. Set by default to \code{NULL}.
#'@param node character string. Node of which the conditional probability distribution is being changed.
#'@param value_node character string. Level of \code{node}.
#'@param value_parents character string. Levels of \code{node}'s parents. The levels should be defined according to the order of the parents in \code{bnfit[[node]][["parents"]]}. If \code{node} has no parents, then should be set to \code{NULL}.
#'@param new_value numeric vector with elements between 0 and 1. Values to which the parameter should be updated. It can take a specific value or more than one. For more than one value, these should be defined through a vector with an increasing order of the elements. \code{new_value} can also take as value the character string \code{all}: in this case a sequence of possible parameter changes ranging from 0.05 to 0.95 is considered.
#'@param covariation character string. Co-variation scheme to be used for the updated Bayesian network. Can take values \code{uniform}, \code{proportional}, \code{orderp}, \code{all}. If equal to \code{all}, uniform, proportional and order-preserving co-variation schemes are considered. Set by default to \code{proportional}.
#'
#'@references Coupé, V. M., & Van Der Gaag, L. C. (2002). Properties of sensitivity analysis of Bayesian belief networks. Annals of Mathematics and Artificial Intelligence, 36(4), 323-356.
#'@references Leonelli, M., Goergen, C., & Smith, J. Q. (2017). Sensitivity analysis in multilinear probabilistic models. Information Sciences, 411, 84-97.
#'
#' @examples sensitivity(synthetic_bn, "y2", "3", node = "y1",value_node = "1",
#'  value_parents = NULL, new_value = "all", covariation = "all")
#' @examples sensitivity(synthetic_bn, "y3", "1", "y2", "1", node = "y1", "1", NULL, 0.9, "all")
#'
#'@import bnlearn
#'@importClassesFrom bnlearn bn.fit
#'@importFrom stats coef
#'@importFrom bnlearn as.grain
#'@importFrom gRain setEvidence
#'@importFrom gRain querygrain
#'@importFrom tidyr gather
#'@importFrom gRbase compile
#'@export
sensitivity <- function(bnfit,
                        interest_node,
                        interest_node_value,
                        evidence_nodes = NULL,
                        evidence_states = NULL,
                        node,
                        value_node,
                        value_parents,
                        new_value,
                        covariation = "proportional"
                        )
{
  if (!(interest_node %in% names(bnfit))) {
    stop("Invalid input for interest_node")
  }
  if(!is.null(evidence_nodes)){
    if (any(!(evidence_nodes %in% names(bnfit)))) {
      stop("Invalid input for evidence_nodes")
    }
  }
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
        "Order preserving co-variation scheme can't work by varying the last parameter in the increasing order of the parameters."
      )
    }
    if (s == length.node & covariation == "all") {
      warning(
        "The sensitivity function with order-preserving scheme won't be returned since the last parameter in the increasing order of the parameters is being varied."
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
      sens <-
        data.frame(new_value2, rep(NA, length(new_value2)), rep(NA, length(new_value2)))
      colnames(sens) <- c("New_value", "Uniform", "Proportional")
    } else{
      sens <-
        data.frame(new_value2,
                   rep(NA, length(new_value2)),
                   rep(NA, length(new_value2)),
                   rep(NA, length(new_value2)))
      colnames(sens) <-
        c("New_value",
          "Uniform",
          "Proportional",
          "Order Preserving")
    }
  } else{
    sens <- data.frame(new_value2, rep(NA, length(new_value2)))
    colnames(sens) <- c("New_value", "Sensitivity")
  }

  if (length(new_value2) == length(new_value_op)) {
    for (k in 1:length(bnfit.new.scheme)) {
      for (j in 1:length(new_value2)) {
     #   suppressWarnings(if (!is.na(bnfit.new.scheme[[k]][[j]])) {
          grain.bn <- gRbase::compile(bnlearn::as.grain(bnfit.new.scheme[[k]][[j]]))
          suppressWarnings(if (!is.null(evidence_nodes)) {
            grain.bn <-
              gRain::setEvidence(grain.bn, nodes = evidence_nodes, states = evidence_states)
          })
          sens[j, k + 1] <-
            gRain::querygrain(grain.bn, nodes = interest_node)[[interest_node]][interest_node_value]
 #       }
  #      )
      }
    }
  } else{
    if (covariation == "all") {
      for (k in 1:(length(bnfit.new.scheme) - 1)) {
        for (j in 1:length(new_value2)) {
          suppressWarnings(if (all(!is.na(bnfit.new.scheme[[k]][[j]]))) {
            grain.bn <- gRbase::compile(bnlearn::as.grain(bnfit.new.scheme[[k]][[j]]))
            if (!is.null(evidence_nodes)) {
              grain.bn <-
                gRain::setEvidence(grain.bn, nodes = evidence_nodes, states = evidence_states)
            }
            sens[j, k + 1] <-
              gRain::querygrain(grain.bn, nodes = interest_node)[[interest_node]][interest_node_value]
          })
        }
      }
      if (length(new_value_op) != 0) {
        for (t in 1:(length(new_value_op))) {
          suppressWarnings(if (all(!is.na(bnfit.new.scheme[[3]][[t]]))) {
            grain.bn <- gRbase::compile(bnlearn::as.grain(bnfit.new.scheme[[3]][[t]]))
            if (!is.null(evidence_nodes)) {
              grain.bn <-
                gRain::setEvidence(grain.bn, nodes = evidence_nodes, states = evidence_states)
            }
            sens[t, 4] <-
              gRain::querygrain(grain.bn, nodes = interest_node)[[interest_node]][interest_node_value]
          })
        }
      }
    } else{
      if (length(new_value_op) != 0) {
        for (t in 1:(length(new_value_op))) {
          suppressWarnings(if (all(!is.na(bnfit.new.scheme[[1]][[t]]))) {
            grain.bn <- gRbase::compile(bnlearn::as.grain(bnfit.new.scheme[[1]][[t]]))
            if (!is.null(evidence_nodes)) {
              grain.bn <-
                gRain::setEvidence(grain.bn, nodes = evidence_nodes, states = evidence_states)
            }
            sens[t, 2] <-
              gRain::querygrain(grain.bn, nodes = interest_node)[[interest_node]][interest_node_value]
          })
        }
      }
    }
  }
 # if (all(is.na(sens[, -1]))) {
#    plot <- FALSE
#    warning("The plot won't be showed since all the values are not possible")
#  }
    if (covariation == "all") {
      if (nrow(sens) == 1) {
        ci <- tidyr::gather(sens, key = "scheme", value = "value", - New_value)
        New_value <- ci$New_value
        value <- ci$value
        scheme <- ci$scheme
        Sensitivity <- ci$Sensitivity
        plot <- ggplot(ci, aes(x = New_value, y = value)) + geom_point(aes(color = scheme)) + ylab("Sensitivity") + theme_minimal() + xlab("New value")
      } else{
        ci <- tidyr::gather(sens, key = "scheme", value = "value", - New_value)
        New_value <- ci$New_value
        value <- ci$value
        scheme <- ci$scheme
        Sensitivity <- ci$Sensitivity
        plot <- ggplot(ci, aes(x = New_value, y = value)) + geom_line(aes(color = scheme)) + ylab("Sensitivity") + theme_minimal() + xlab("New value")
      }
    } else{
      if (nrow(sens) == 1) {
        plot <- ggplot(sens, aes(x = New_value, y = Sensitivity)) +geom_point() + theme_minimal()
      } else{
        plot <- ggplot(sens,aes(x = New_value, y = Sensitivity)) + geom_line() + theme_minimal()
      }
    }
  out <- list(sensitivity = sens, plot = plot)
  attr(out,'class') <- 'sensitivity'
  return(out)
}



#@examples sensitivity(synthetic_bn, "y2", "3", node = "y1",value_node = "1", value_parents = NULL,
#    new_value = "all", covariation = "all")
#@examples sensitivity(synthetic_bn, "y3", "1", "y2", "1", node = "y1", "1", NULL, 0.9, "all")


