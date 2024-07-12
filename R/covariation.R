#' Co-variation schemes
#'
#' Functions that return an updated Bayesian network using the proportional, uniform and order-preserving co-variation schemes.
#'
#' The Bayesian network on which parameter variation is being conducted should be expressed as a \code{bn.fit} object.
#' The name of the node to be varied, its level and its parent's levels should be specified.
#' The parameter variation specified by the function is:
#'
#'  P ( \code{node} = \code{value_node} | parents = \code{value_parents} ) = \code{new_value}
#'
#' For \code{orderp_covar}, if two or more parameters in a distribution have the same value, the order is given by the one in the respective conditional probability table. Furthermore, the parameter associated to the largest probability of the conditional probability law cannot be varied.
#'
#' @return An object of class \code{bn.fit} with updated probabilities.
#'
#'@param bnfit object of class \code{bn.fit}.
#'@param node character string. Node of which the conditional probability distribution is being changed.
#'@param value_node character string. Level of \code{node}.
#'@param value_parents character string. Levels of \code{node}'s parents. The levels should be defined according to the order of the parents in \code{bnfit[[node]][["parents"]]}. If \code{node} has no parents, then it should be set to \code{NULL}.
#'@param new_value numeric value between 0 and 1. Value to which the parameter should be updated.
#'
#'@references Laskey, K. B. (1995). Sensitivity analysis for probability assessments in Bayesian networks. IEEE Transactions on Systems, Man, and Cybernetics, 25(6), 901-909.
#'@references Renooij, S. (2014). Co-variation for sensitivity analysis in Bayesian networks: Properties, consequences and alternatives. International Journal of Approximate Reasoning, 55(4), 1022-1042.
#'@references Leonelli, M., & Riccomagno, E. (2022). A geometric characterization of sensitivity analysis in monomial models. International Journal of Approximate Reasoning, 151, 64-84.
#'#'
#'@examples proportional_covar(synthetic_bn, "y3", "2", c("2","1"), 0.3)
#'@examples uniform_covar(synthetic_bn, "y2", "1", "2", 0.3)
#'@examples orderp_covar(synthetic_bn, "y1", "1", NULL, 0.3)
#'
#'@name covariation
#'@importFrom stats coef
#'@importClassesFrom bnlearn bn.fit
NULL

#' @rdname covariation
#' @importFrom stats coef
#'@importClassesFrom bnlearn bn.fit
#'@export
#'
proportional_covar <-
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
    for (i in 1:length.node) {
      b <-
        append(unlist(dimnames(coef(bnfit[[node]]))[1])[i], value_parents)
      cpt.new[t(b)] <- (1 - new_value) / (1 - old) * cpt.new[t(b)]
    }
    cpt.new[t(append(value_node, value_parents))] <- new_value
    bnfit[[node]] <- cpt.new
    return(bnfit)
  }

#' @rdname covariation
#' @importFrom stats coef
#'@importClassesFrom bnlearn bn.fit
#'@export
#'
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
        "Order preserving co-variation scheme can't work by varying the last parameter in the increasing order of the parameters"
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

#' @rdname covariation
#' @importFrom stats coef
#'@importClassesFrom bnlearn bn.fit
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

