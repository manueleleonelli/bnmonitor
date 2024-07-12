#' Sensitivity of probability query
#'
#' \code{sensquery} returns, for a given change in a probability of interest, the parameters' changes to achieve it together with the corresponding CD distances.
#'
#'
#' The Bayesian network should be expressed as a \code{bn.fit} object.
#' The name of the node of the probability of interest, its level and the new value should be specified. Evidence could be also indicated.
#' The probability of interest is specified as follows:
#'
#'  P ( \code{interest_node} = \code{interest_node_value} | \code{evidence_nodes} = \code{evidence_states}  ) = \code{new_value}
#'
#' Only the  proportional co-variation scheme is used.
#'
#' @seealso \code{\link{sensitivity}}
#'
#' @return A dataframe with the following columns: \code{node} - the vertex of the proposed change; \code{Value node} - the level of \code{node} to be changed; \code{Value parents} - the levels of the parent variables of \code{node}; \code{Original value} - the original probability defined by \code{Node}, \code{Value node} and \code{Value parents}; \code{Suggested change} - the new proposed value for the probability defined by \code{Node}, \code{Value node} and \code{Value parents}; \code{CD distance} - the CD distance between the original and new network with the \code{Suggested change}.
#'
#'@param bnfit object of class \code{bn.fit}.
#'@param interest_node character string. Node of the probability query of interest.
#'@param interest_node_value character string. Level of \code{interest_node}.
#'@param new_value numeric value between 0 and 1. New value of the probability of interest.
#'@param evidence_nodes character string. Evidence nodes. Set by default to \code{NULL}.
#'@param evidence_states character string. Levels of \code{evidence_nodes}. If \code{NULL} no evidence is considered. If \code{evidence_nodes="NULL"}, \code{evidence_states} should be set to \code{NULL}. Set by default to \code{NULL}.
#'
#'@examples sensquery(synthetic_bn,"y3", "3", 0.3)
#'
#'@references Chan, H., & Darwiche, A. (2002). When do numbers really matter?. Journal of Artificial Intelligence Research, 17, 265-287.
#'
#'@import bnlearn
#'@importClassesFrom bnlearn bn.fit
#'@importFrom stats coef complete.cases uniroot
#'@export
sensquery <- function(bnfit,
                      interest_node,
                      interest_node_value,
                      new_value,
                      evidence_nodes = NULL,
                      evidence_states = NULL) {

  if (!(interest_node %in% names(bnfit))) {
    stop("Invalid input for interest_node")
  }
  if(!is.null(evidence_nodes)){
    if (any(!(evidence_nodes %in% names(bnfit)))) {
      stop("Invalid input for evidence_nodes")
    }
  }
  query <- data.frame()
  for (i in 1:length(bnfit)) {
    original <-
      numeric(nrow(expand.grid(dimnames(bnfit[[i]][["prob"]]))))
    i_change <-
      numeric(nrow(expand.grid(dimnames(bnfit[[i]][["prob"]]))))
    ij_parents <-
      character(nrow(expand.grid(dimnames(bnfit[[i]][["prob"]]))))
    for (j in 1:nrow(expand.grid(dimnames(bnfit[[i]][["prob"]])))) {
      value_node <-
        as.character.factor(expand.grid(dimnames(bnfit[[i]][["prob"]]))[j, 1])
      if (length(bnfit[[i]][["parents"]]) == 0) {
        value_parents <- NULL
      } else{
        value_parents <-
          as.character.factor(unname(unlist(expand.grid(
            dimnames(bnfit[[i]][["prob"]])
          )[j, -1])))
      }
      ij_parents[j] <- paste(value_parents, collapse = ",")
      original[j] <-
        coef(bnfit[[i]])[t(append(value_node, value_parents))]
      i_change[j] <- try(uniroot(
        fquery,
        c(0, 1),
        new_value = new_value,
        bnfit = bnfit,
        interest_node = interest_node,
        interest_node_value = interest_node_value,
        evidence_nodes = evidence_nodes,
        evidence_states = evidence_states,
        node = nodes(bnfit)[i],
        value_node = value_node,
        value_parents = value_parents,
        extendInt = "yes"
      )$root,
      silent = TRUE)
    }
    query <-
      suppressWarnings(rbind(
        query,
        data.frame(
          node = bnfit[[i]]$node,
          value_node = expand.grid(dimnames(bnfit[[i]][["prob"]]))[, 1],
          value_parents = ij_parents,
          original = as.numeric(original),
          change = as.numeric(i_change)
        )
      ))
  }
  final_query <- query[complete.cases(query), ]
  if (nrow(final_query) == 0) {
    final_query <- NULL
    warning("The BN satisfies the given constraint. No parameter changes are needed.")
  } else{
    final_query <- subset(final_query, !duplicated(final_query[, 1]))
    dist <- numeric(nrow(final_query))
    for (t in 1:nrow(final_query)) {
      dist[t] <-
        CD(
          bnfit = bnfit,
          node = as.character(final_query[t, 1]) ,
          value_node = as.character(final_query[t, 2]),
          value_parents = unlist(strsplit(as.character(final_query[t, 3]), ",")),
          new_value = final_query[t, 5],
          covariation = "proportional"
        )$CD[, 2]
    }
    final_query[, 6] <- dist
    final_query <- final_query[order(final_query[, 6]), ]
    rownames(final_query) <- seq(1, nrow(final_query))
    colnames(final_query) <-
      c(
        "Node",
        "Value node",
        "Value parents",
        "Original value",
        "Suggested change",
        "CD distance"
      )
  }
  return(final_query)
}
