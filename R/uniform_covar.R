#' Uniform co-variation scheme
#'
#' \code{uniform_covar} returns an updated Bayesian network using uniform co-variation scheme.
#'
#' \cr
#' The Bayesian network on which parameter variation is being conducted should be expressed as a bn.fit object.
#' The name of the node to be varied, its level and its parent's level should be specified.
#' The parameter variation specified by the function is:\cr
#' \cr
#' P ( \code{node} = \code{value_node} | parents = \code{value_parents} ) = \code{new_value} \cr
#'
#' @family co-variation schemes
#'
#'@param bnfit object of class bn.fit
#'@param node character string. Node of which the conditional probability distribution is being changed.
#'@param value_node character string. Level of \code{node}.
#'@param value_parents character string. Levels of \code{node}'s parents. The levels should be defined according to the order of the parents in \code{bnfit[[node]][["parents"]]}. If \code{node} has no parents, then should be set to \code{NULL}.
#'@param new_value numeric value between 0 and 1. Value to which the parameter should be updated.
#'
#'@importClassesFrom bnlearn bn.fit
#'@importFrom stats coef
#'@export
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
