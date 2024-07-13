#' Strength of edges in a Bayesian network
#'
#' Computation of the measure of edge strength for all edges in a Bayesian networks
#'
#' The measure of edge strength is defined as the largest diameter out of all conditional probability tables where all other parents but the considered one are fixed to a specific combination.
#'
#' @return A dataframe with first two columns the edge list of the \code{bn.fit} input object. The third column \code{Edge.Strength} reports the measure of edge strength.
#'
#'@param bnfit object of class \code{bn.fit}.
#'
#' @seealso \code{\link{diameter}}
#'
#'@examples edge_strength(travel)
#'
#'@references Leonelli, M., Smith, J. Q., & Wright, S. K. (2024). The diameter of a stochastic matrix: A new measure for sensitivity analysis in Bayesian networks. arXiv preprint arXiv:2407.04667.
#'
#'@import bnlearn
#'@importClassesFrom bnlearn bn.fit
#'@importFrom bnlearn nodes root.nodes arcs
#'@export


edge_strength <- function(bnfit){
  if(class(bnfit)[1] != "bn.fit"){stop("The input is not a bn.fit object")}
  edge_strength <- rep(0,nrow(arcs(bnfit)))
  for(m in 1:nrow(arcs(bnfit))){
    object <- bnfit[[unname(arcs(bnfit)[m,2])]]
    s <- which(object$parents == unname(arcs(bnfit)[m,1]))
    if(s > 1){
      entry <- s
      tot_parents <- length(object$parents)
      levels_parents <- unname(sapply(dimnames(object$prob),function(i) length(i)))[-1]
      lev <- levels_parents[entry]
      lev_node <- unname(sapply(dimnames(object$prob),function(i) length(i)))[1]
      number <- lev*lev_node
      breaks <- prod(levels_parents[entry:tot_parents])/levels_parents[entry]
      tot <- length(object$prob)
      modulo <- lev_node*prod(levels_parents[1:entry])/levels_parents[entry]
      indexes <- c(which((1:tot%%(tot/breaks)) == 1),tot+1)
      probi <- rep(0,tot)
      times <- (indexes[2]-indexes[1])/number
      swap <- c()
      for(i in 1:(length(indexes)-1)){
        original <- object$prob[indexes[i]:indexes[i+1]]
        for(j in 1:times){
          for(k in 1:lev){
            swap <- c(swap,((((j-1)*lev_node + 1)+((k-1)*modulo))+indexes[i]-1)    : ((((j-1)*lev_node + 1)+((k-1)*modulo ) + lev_node - 1)+indexes[i]-1))
          }
        }
      }
      probi <- object$prob[swap]
      dim(probi) <- c(dim(object$prob)[1],dim(object$prob)[entry+1],dim(object$prob)[-c(1,entry+1)])
      parents <- c(object$parents[entry],object$parents[-entry])
    }
    if(s == 1){
      parents <- object$parents
      probi <- object$prob
    }
    levels_out <- dim(probi)[1]
    levels_parent <- dim(probi)[2]
    tot <- length(probi)
    indexes <- c(which(1:tot%%(levels_out*levels_parent) == 1),tot)
    ind <- which(1:tot%% levels_out == 1)
    max <- 0
    for(i in 1:(length(indexes)-1)){
      values <- ind[ind>= indexes[i] & ind < indexes[i+1]]
      for(j in 1:(length(values)-1)){
        for(k in (j+1):length(values)){
          temp <- tvd(probi[values[j]:(values[j]+levels_out-1)],probi[values[k]:(values[k]+levels_out-1)])
          if(max < temp) max <- temp
        }
      }
    }
    edge_strength[m] <- max
  }
  result <- data.frame(arcs(bnfit),Edge_Strength = edge_strength)
  attr(result, 'class') <- c('edgestrength','data.frame')
  return(result)
}

