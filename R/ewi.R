#' Edge-weigthed influence
#'
#' @description
#' Computation of the edge-weigthed influence in a Bayesian network
#'
#' @details
#' The edge-weigthed influence of a node \eqn{X_j} on an output node \eqn{X_i} in a Bayesian network is \deqn{EWI(X_j,X_i)= \sum_{s\in S_{ji}}\left(\prod_{(k,l)\in s}\delta_{kl}\right)^{|s|},} where \eqn{S_{ji}} is the set of active trails between \eqn{X_j} and \eqn{X_i}, \eqn{\delta_{kl}} is the strength of an edge between \eqn{X_k} and \eqn{X_l}, and \eqn{|s|} is the length of the trail \eqn{s}.
#'
#'
#' @return A dataframe with the following columns: \code{Nodes} - the vertices of the BN; \code{Influence} - the edge-weigthed influence of the corresponding node.
#'
#'@param bnfit object of class \code{bn.fit}.
#'@param node a node of \code{bnfit}
#'
#'@examples ewi(travel, "T")
#'
#'@seealso \code{\link{mutual_info}}, \code{\link{dwi}}, \code{\link{edge_strength}}
#'
#'@references Leonelli, M., Smith, J. Q., & Wright, S. K. (2024). The diameter of a stochastic matrix: A new measure for sensitivity analysis in Bayesian networks. arXiv preprint arXiv:2407.04667.
#'
#'@import bnlearn
#'@importClassesFrom bnlearn bn.fit
#'@importFrom bnlearn nodes root.nodes parents as.igraph
#'@importFrom igraph all_simple_paths
#'@export

ewi <- function(bnfit, node){
  if(class(bnfit)[1] != "bn.fit"){stop("An input bn.fit object required")}
  if(!node%in% nodes(bnfit)){stop("Not a valid vertex name in input")}
  pos.nodes <- nodes(bnfit)[-which(nodes(bnfit)== node)]
  strength <- edge_strength(bnfit)
  weigth <- rep(0,length(pos.nodes))
  dag_temp <- as.igraph(bnfit)
  for(k in 1:length(pos.nodes)){
    ciao <- lapply(all_simple_paths(dag_temp,node,pos.nodes[k],mode="all"),names)
    if(length(ciao) != 0){
      for(i in 1:length(ciao)){
        if(length(ciao[[i]]) == 2) {weigth[k] <-  weigth[k] + strength[which( (strength[,1] == ciao[[i]][1] & strength[,2] == ciao[[i]][2]) | (strength[,1] == ciao[[i]][2] & strength[,2] == ciao[[i]][1])),3]} else{
          status <- T
          for(j in 3:length(ciao[[i]])){
            if(ciao[[i]][j-2]%in%bnlearn::parents(bnfit,ciao[[i]][j-1]) & ciao[[i]][j]%in%bnlearn::parents(bnfit,ciao[[i]][j-1])) {status <- F}
          }
          if(status){
            temp <- 1
            for(j in 2:length(ciao[[i]])){
              temp <- temp*strength[which( (strength[,1] == ciao[[i]][j-1] & strength[,2] == ciao[[i]][j]) | (strength[,1] == ciao[[i]][j] & strength[,2] == ciao[[i]][j-1])),3]
            }
            weigth[k] <- weigth[k] + temp^(length(ciao[[i]])-1)
          }
        }
      }
    }
  }
  result <- data.frame(Nodes = append(pos.nodes, node, after=which(nodes(bnfit)== node)-1), Influence = append(weigth,NA,after=which(nodes(bnfit)== node)-1))
  result <- list(EWI = result, BN = bnfit)
  attr(result, 'class') <- 'ewi'
  return(result)
  }


