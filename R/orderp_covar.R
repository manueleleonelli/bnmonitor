#' Order-preserving covariation scheme
#'
#' \code{orderp_covar} returns an updated Bayesian network using order-preserving covariation scheme
#'
#' The Bayesian network on which parameter variation is being conducted should be expressed as a bn.fit object.
#' The name of the node to be varied, its level and its parent's level should be specified.
#' The parameter variation specified by the function is:
#'
#'  P ( \code{node} = \code{value_node} | parents = \code{value_parents} ) = \code{new_value}
#'
#'It should be noted that if two or more parameters in a distribution have the same value, the order is given by the one in the respective conditional probability table.
#'
#'@family covariation schemes
#'
#'@param bnfit object of class bn.fit
#'@param node character string. Node of which the conditional probability distribution is being changed.
#'@param value_node character string. Level of \code{node}.
#'@param value_parents character string. Levels of \code{node}'s parents. The levels should be defined according to the order of the parents in \code{bnfit[[node]][["parents"]]}. If \code{node} has no parents, then should be set to \code{NULL}.
#'@param new_value numeric value between 0 and 1. Value to which the parameter should be updated.
#'

#'@importFrom stats coef
#'@importClassesFrom bnlearn bn.fit
#'@export
orderp_covar <-
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
    length.node <- dim(coef(bnfit[[node]]))[1]
    old <- coef(bnfit[[node]])[t(append(value_node, value_parents))]
    cpt.new <- coef(bnfit[[node]])
    node_values <- numeric(length.node)
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
    if (s == length.node) {
      stop (
        "Order preserving covariation scheme can't work by varying the last parameter in the increasing order of the parameters"
      )
    }
    if (new_value > 1 / (1 + dim(coef(bnfit[[node]]))[1] - s)) {
      bnfit <- NA
      return(bnfit)
    }
    m.T <-
      sum(node_values[order > s])
    x.T <- 1 / (1 + length.node - s)
    for (i in 1:length.node) {
      k <- order[i]
      if (k == s) {
        cpt.new[t(append(unlist(dimnames(coef(
          bnfit[[node]]
        ))[1])[i], value_parents))] <- new_value
      }
      if (k < s) {
        if (new_value <= old) {
          cpt.new[t(append(unlist(dimnames(
            coef(bnfit[[node]])
          )[1])[i], value_parents))] <-
            cpt.new[t(append(unlist(dimnames(
              coef(bnfit[[node]])
            )[1])[i], value_parents))] / old * new_value
        } else{
          cpt.new[t(append(unlist(dimnames(
            coef(bnfit[[node]])
          )[1])[i], value_parents))] <-
            cpt.new[t(append(unlist(dimnames(
              coef(bnfit[[node]])
            )[1])[i], value_parents))] / (x.T - old) * (x.T - new_value)
        }
      }
      if (k > s) {
        if (new_value <= old) {
          cpt.new[t(append(unlist(dimnames(
            coef(bnfit[[node]])
          )[1])[i], value_parents))] <-
            (-cpt.new[t(append(unlist(dimnames(
              coef(bnfit[[node]])
            )[1])[i], value_parents))] * (1 - m.T)) / (m.T * old) * new_value + cpt.new[t(append(unlist(dimnames(
              coef(bnfit[[node]])
            )[1])[i], value_parents))] / m.T
        } else{
          cpt.new[t(append(unlist(dimnames(
            coef(bnfit[[node]])
          )[1])[i], value_parents))] <-
            (cpt.new[t(append(unlist(dimnames(
              coef(bnfit[[node]])
            )[1])[i], value_parents))] - x.T) / (x.T - old) * (x.T - new_value) + x.T
        }
      }
    }
    bnfit[[node]] <- cpt.new
    return(bnfit)
  }
