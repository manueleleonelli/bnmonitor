#' Distance-weigthed influence
#'
#' @description
#' Computation of the distance-weigthed influence in a Bayesian network
#'
#' @details
#' The distance-weigthed influence of a node \eqn{X_j} on an output node \eqn{X_i} in a Bayesian network is \deqn{DWI(X_j,X_i,w)= \sum_{s\in S_{ji}}w^{|s|},} where \eqn{S_{ji}} is the set of active trails between \eqn{X_j} and \eqn{X_i}, \eqn{w\in(0,1]} is an input parameter, and \eqn{|s|} is the length of the trail \eqn{s}.
#'
#'
#' @return A dataframe with the following columns: \code{Nodes} - the vertices of the BN; \code{Influence} - the distance-weigthed influence of the corresponding node.
#'
#'@param bn object of class \code{bn.fit} or \code{bn}.
#'@param node a node of \code{bnfit}.
#'@param w a number in \eqn{(0,1]}.
#'
#'@seealso \code{\link{ewi}}, \code{\link{mutual_info}}
#'
#'@examples dwi(travel, "T", 0.5)
#'
#'@references Albrecht, D., Nicholson, A. E., & Whittle, C. (2014). Structural sensitivity for the knowledge engineering of Bayesian networks. In Probabilistic Graphical Models (pp. 1-16). Springer International Publishing.
#'
#'@import bnlearn
#'@importClassesFrom bnlearn bn.fit
#'@importFrom bnlearn nodes root.nodes parents as.igraph bn.net
#'@importFrom igraph all_simple_paths
#'@importFrom methods is
#'@export

dwi <- function(bn, node, w){
  if(class(bn)[1] == "bn.fit"){bn <- bn.net(bn)}
  if(!is(bn,"bn")){stop("Not a valid bnlearn object in input")}
  if(!node%in% nodes(bn)){stop("Not a valid vertex name in input")}
  if(w>1 | w <=0){stop("The input w must be between zero and one")}
  pos.nodes <- nodes(bn)[-which(nodes(bn)== node)]
  weigth <- rep(0,length(pos.nodes))
  dag_temp <- as.igraph(bn)
  for(k in 1:length(pos.nodes)){
    ciao <- lapply(all_simple_paths(dag_temp,node,pos.nodes[k],mode="all"),names)
    if(length(ciao) != 0){
      for(i in 1:length(ciao)){
        if(length(ciao[[i]]) == 2) {weigth[k] <-  weigth[k] + w} else{
          status <- T
          for(j in 3:length(ciao[[i]])){
            if(ciao[[i]][j-2]%in%bnlearn::parents(bn,ciao[[i]][j-1]) & ciao[[i]][j]%in%bnlearn::parents(bn,ciao[[i]][j-1])) {status <- F}
          }
          if(status){weigth[k] <- weigth[k] + w^(length(ciao[[i]])-1)}
        }
      }
    }
  }
  result <- data.frame(Nodes = append(pos.nodes, node, after=which(nodes(bn)== node)-1), Influence = append(weigth,NA,after=which(nodes(bn)== node)-1))
  result <- list(DWI = result, BN = bn)
  attr(result, 'class') <- 'dwi'
  return(result)
}

