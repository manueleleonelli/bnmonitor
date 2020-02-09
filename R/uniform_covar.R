#' Uniform covariation scheme
#'
#' \code{uniform_covar} returns an updated Bayesian network using the uniform covariation scheme.
#'
#' The Bayesian network on which parameter variation is being conducted should be expressed as a \code{bn.fit} object.
#' The name of the node to be varied, its level and its parent's levels should be specified.
#' The parameter variation specified by the function is:
#'
#' P ( \code{node} = \code{value_node} | parents = \code{value_parents} ) = \code{new_value}
#'
#'@family covariation schemes
#'
#'@param bnfit object of class \code{bn.fit}.
#'@param node character string. Node of which the conditional probability distribution is being changed.
#'@param value_node character string. Level of \code{node}.
#'@param value_parents character string. Levels of \code{node}'s parents. The levels should be defined according to the order of the parents in \code{bnfit[[node]][["parents"]]}. If \code{node} has no parents, then it should be set to \code{NULL}.
#'@param new_value numeric value between 0 and 1. Value to which the parameter should be updated.
#'
#'@references Renooij, S. (2014). Co-variation for sensitivity analysis in Bayesian networks: Properties, consequences and alternatives. International journal of approximate reasoning, 55(4), 1022-1042.
#'@references Leonelli, M., & Riccomagno, E. (2018). A geometric characterisation of sensitivity analysis in monomial models. arXiv preprint arXiv:1901.02058.
#'
#'@examples uniform_covar(synthetic_bn, "y3", "2", c("2","1"), 0.3)
#'@examples uniform_covar(synthetic_bn, "y2", "1", "2", 0.3)
#'@examples uniform_covar(synthetic_bn, "y1", "1", NULL, 0.3)
#'
#'@importClassesFrom bnlearn bn.fit
#'@importFrom stats coef
#'@export
#'
uniform_covar <-
  function(bnfit,
           node,
           value_node,
           value_parents,
           new_value) {
    if (!(node %in% names(bnfit))) {
      stop("Invalid input for node")
    }
    nparents <- length(bnfit[[node]][["parents"]])
    if (length(value_parents) != nparents) {
      stop("Invalid length of value_parents")
    }
    cpt.new <- coef(bnfit[[node]])
    length.node <- dim(coef(bnfit[[node]]))[1]
    for (i in 1:length.node) {
      b <- append(unlist(dimnames(coef(bnfit[[node]]))[1])[i], value_parents)
      cpt.new[t(b)] <- (1 - new_value) / (length.node - 1)
    }
    cpt.new[t(append(value_node, value_parents))] <- new_value
    bnfit[[node]] <- cpt.new
    return(bnfit)
  }
