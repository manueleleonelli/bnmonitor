#' Mutual information
#'
#' @description
#' Computation of the mutual information in a Bayesian network
#'
#' @details
#' The mutual information between two variables \eqn{X_j} and \eqn{X_i} with sample spaces \eqn{\mathcal{X}_i} and \eqn{\mathcal{X}_j}, respectively, is equal to \deqn{\sum_{x_j\in\mathcal{X}_j}\sum_{x_i\in\mathcal{X}_i}p(x_i,x_j)\log\frac{p(x_i,x_j)}{p(x_i)p(x_j)}.}
#'
#' @return A dataframe with the following columns: \code{Nodes} - the vertices of the BN; \code{MutualInfo} - the mutual information of the corresponding node.
#'
#' @seealso \code{\link{ewi}}, \code{\link{dwi}}
#'
#'@param bnfit object of class \code{bn.fit}.
#'@param node a node of \code{bnfit}.
#'
#'@examples mutual_info(travel, "T")
#'
#'@references Albrecht, D., Nicholson, A. E., & Whittle, C. (2014). Structural sensitivity for the knowledge engineering of Bayesian networks. In Probabilistic Graphical Models (pp. 1-16). Springer International Publishing.
#'
#'@import bnlearn
#'@importClassesFrom bnlearn bn.fit
#'@importFrom bnlearn nodes root.nodes parents as.grain bn.net
#'@importFrom gRain querygrain
#'@export


mutual_info <- function(bnfit, node){
  if(class(bnfit)[1] != "bn.fit"){stop("An input bn.fit object required")}
  if(!node%in% nodes(bnfit)){stop("Not a valid vertex name in input")}
  possible <- nodes(bnfit)[-which(nodes(bnfit)==node)]
  result <- data.frame(Nodes = append(possible, node, after=which(nodes(bnfit)== node)-1), MutalInfo = append(unname(sapply(possible,mi, bnfit= bnfit,node1 = node)),NA,after=which(nodes(bnfit)== node)-1)  )
  result <- list(MutualInfo = result, BN = bnfit)
  attr(result, 'class') <- 'mutualinfo'
  return(result)
}

