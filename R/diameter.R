#' Diameters in a Bayesian network
#'
#' Computation of the diameters of all conditional probability tables in a Bayesian network.
#'
#'
#' The diameter of a conditional probability table \eqn{P} with \eqn{n} rows \eqn{p_1,\dots,p_n} is \deqn{d^+(P)=\max_{i,j\leq n} d_V(p_i,p_j),} where \eqn{d_V} is the total variation distance between two probability mass functions over a sample space \eqn{\mathcal{X}}, i.e. \deqn{d_V(p_i,p_j)=\frac{1}{2}\sum_{x\in\mathcal{X}}|p_i(x)-p_j(x)|.}
#'
#'
#' @return A dataframe with the following columns: \code{Nodes} - the vertices of the BN; \code{Diameter} - the diameters of the associated conditional probability tables.
#'
#'@param bnfit object of class \code{bn.fit}.
#'
#'@examples diameter(travel)
#'
#'@references Leonelli, M., Smith, J. Q., & Wright, S. K. (2024). The diameter of a stochastic matrix: A new measure for sensitivity analysis in Bayesian networks. arXiv preprint arXiv:2407.04667.
#'
#'@import bnlearn
#'@importClassesFrom bnlearn bn.fit
#'@importFrom bnlearn nodes root.nodes
#'@export

diameter<- function(bnfit){
  if(class(bnfit)[1] != "bn.fit"){stop("The input is not a bn.fit object")}
  out <- rep(NA,length(nodes(bnfit)))
  for(i in 1:length(nodes(bnfit))){
    if(!nodes(bnfit)[i]%in% root.nodes(bnfit)){
      out[i] <- diam(bnfit,nodes(bnfit)[i])
    }
  }
  result <- data.frame(Nodes = nodes(bnfit),Diameter = unname(out))
  result <- list(Diameter = result, BN = bnfit)
  attr(result, 'class') <- 'diameter'
  return(result)
}

