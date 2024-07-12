#' Measures of asymmetric independence
#'
#' @description
#'  Computation of the indexes of context-specific and partial independence in a conditional probability table.
#'
#' @details
#' The index of context-specific independence computes the upper diameter of a CPT where all parents but one are fixed, while the index of partial independence computes the lower diameter for the same CPT. If the lower diameter is close to zero it means that there are at least two rows of the CPT which are very similar to each other, thus implying a partial conditional independence.
#'
#'
#' @return A list where each entry refers to one of the parent of \code{node}. Each entry of the list is a dataframe with the index of context-specific independence for each outcome of the conditioning parent in one column and in additional column the index of partial independence in the case of non-binary variables.
#'
#' @param bnfit object of class \code{bn.fit}.
#' @param node a node of \code{bnfit}.
#'
#' @seealso \code{\link{diameter}}
#'
#'@examples asy_measure(synthetic_bn, "y3")
#'
#'@references Leonelli, M., Smith, J. Q., & Wright, S. K. (2024). The diameter of a stochastic matrix: A new measure for sensitivity analysis in Bayesian networks. arXiv preprint arXiv:2407.04667.
#'@references Pensar, J., Nyman, H., Lintusaari, J., & Corander, J. (2016). The role of local partial independence in learning of Bayesian networks. International Journal of Approximate Reasoning, 69, 91-105.
#'
#'@import bnlearn
#'@importClassesFrom bnlearn bn.fit
#'@importFrom bnlearn nodes root.nodes arcs parents
#'@export


asy_measure <- function(bnfit, node){
  if(class(bnfit)[1] != "bn.fit"){stop("The input is not a bn.fit object")}
  if(!node%in% nodes(bnfit)){stop("Not a valid vertex name in input")}
  ps <- bnlearn::parents(bnfit,node)
  if(length(ps) < 2){stop("Asymmetry measures only valid for nodes with at least two parents")}
  ls <-  vector("list", length(ps))
  object <- bnfit[[node]]
  for(s in 1:length(ps)){
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
      dimnames(probi)[[1]] <- dimnames(object$prob)[[1]]
      dimnames(probi)[[2]] <- dimnames(object$prob)[[entry+1]]
      pri <- dimnames(object$prob)
      pri[[1]] <- NULL
      pri[[entry]] <- NULL
      if(length(dim(probi))>2){
        for(j in 1:length(pri)) dimnames(probi)[[2+j]] <- pri[[j]]
      }
      parents <- c(object$parents[entry],object$parents[-entry])
      names(dimnames(probi)) <- c(names(dimnames(object$prob))[1],parents)
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
    yy <- dimnames(probi)
    yy[[1]] <- NULL
    yy[[1]] <- NULL
    out <- expand.grid(yy)
    out$Context <- rep(0,nrow(out))
    for(i in 1:(length(indexes)-1)){
      max <- 0
      values <- ind[ind>= indexes[i] & ind < indexes[i+1]]
      for(j in 1:(length(values)-1)){
        for(k in (j+1):length(values)){
          temp <- tvd(probi[values[j]:(values[j]+levels_out-1)],probi[values[k]:(values[k]+levels_out-1)])
          if(max < temp) max <- temp
        }
      }
      out$Context[i] <- max
    }
    if(dim(probi)[2]>2){
      out$Partial <- rep(0,nrow(out))
      for(i in 1:(length(indexes)-1)){
        max <- 100
        values <- ind[ind>= indexes[i] & ind < indexes[i+1]]
        for(j in 1:(length(values)-1)){
          for(k in (j+1):length(values)){
            temp <- tvd(probi[values[j]:(values[j]+levels_out-1)],probi[values[k]:(values[k]+levels_out-1)])
            if(max > temp) max <- temp
          }
        }
        out$Partial[i] <- max
      }
    }
    ls[[s]] <- out
  }
  return(ls)
}
