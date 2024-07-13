#' Amalgamation of levels
#'
#' @description
#' Computation of the diameter of children conditional probability tables when levels of a variable are merged.
#'
#'
#' @return A list where each entry refers to a child of \code{node}. For each child, the function reports the diameter resulting from the amalgamation of any pair of levels.
#'
#'@param bnfit object of class \code{bn.fit}.
#'@param node a node of \code{bnfit}
#'
#'@examples amalgamation(synthetic_bn, "y1")
#'
#'
#'@references Leonelli, M., Smith, J. Q., & Wright, S. K. (2024). The diameter of a stochastic matrix: A new measure for sensitivity analysis in Bayesian networks. arXiv preprint arXiv:2407.04667.
#'
#'@import bnlearn
#'@importClassesFrom bnlearn bn.fit
#'@importFrom bnlearn nodes root.nodes parents children
#'@export

amalgamation <- function(bnfit, node){
  if(class(bnfit)[1] != "bn.fit"){stop("An input bn.fit object required")}
  if(!node%in% nodes(bnfit)){stop("Not a valid vertex name in input")}
  pos.nodes <- bnlearn::children(bnfit,node)
  if(length(pos.nodes) < 1){stop("The chosen node must have children")}
  ls <-  vector("list", length(pos.nodes))
  objects <- dimnames(bnfit[[node]]$prob)[[1]]
  if(length(objects)<3) stop("The chosen node must have at least three levels")
  nm <- rep(NA,choose(length(objects),2))
  k <- 1
  for(i in 1:(length(objects)-1)){
    for (j in (i+1):length(objects)){
      nm[k] <- paste0(objects[i]," - ",objects[j])
      k <- k+1
    }
  }
  for(i in 1:length(pos.nodes)) {
    ls[[i]] <- rep(0,length(nm))
    names(ls[[i]]) <- nm
  }
  for(r in 1:length(pos.nodes)){
    object <- bnfit[[pos.nodes[r]]]
    entrance <- bnlearn::parents(bnfit,pos.nodes[r])
    if(entrance[1] == node){
      parents <- object$parents
      probi <- object$prob
    } else{
      entry <- which(entrance == node)
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
    levels_out <- dim(probi)[1]
    levels_parent <- dim(probi)[2]
    tot <- length(probi)
    indexes <- c(which(1:tot%%(levels_out*levels_parent) == 1),tot)
    ind <- which(1:tot%% levels_out == 1)
    count <- 1
    for(i in 1:(levels_parent-1)){
      for(j in (i+1):levels_parent){
        level1 <- (j + levels_parent*0:length(ind))[(j + levels_parent*0:length(ind)) <= length(ind)]
        level2 <- (i + levels_parent*0:length(ind))[(i + levels_parent*0:length(ind)) <= length(ind)]
        new_probab <- c()
        inde <- indexes
        inde[length(inde)] <- inde[length(inde)]+1
        for(y in 1:(length(inde)-1)){
          new_probab <- c(new_probab,apply(cbind(probi[ind[level1[y]]:(ind[level1[y]]+(levels_out-1))],probi[ind[level2[y]]:(ind[level2[y]]+(levels_out-1))]),1,mean)
                          , probi[inde[y]:(inde[y+1]-1)][!(inde[y]:(inde[y+1]-1)) %in%  c(ind[level1[y]]:(ind[level1[y]]+(levels_out-1)),ind[level2[y]]:(ind[level2[y]]+(levels_out-1))) ]
          )
        }
        dims <- dim(probi)
        dims[2] <- dims[2]-1
        dim(new_probab) <- dims
        names(new_probab) <- names(probi)

        levels_out <- dim(new_probab)[1]
        tot <- length(new_probab)
        indu <- c(which(1:tot%% levels_out == 1),tot+1)
        max <- 0
        for(o in 1:(length(indu)-2)){
          for(p in (1+o):(length(indu)-1)){
            temp <- tvd(new_probab[indu[o]:(indu[o+1]-1)],new_probab[indu[p]:(indu[p+1]-1)])
            if(max < temp) max <- temp
          }
        }
        ls[[r]][count] <- max
        count <- count + 1
      }
    }
  }
  names(ls) <- bnlearn::children(bnfit,node)
  return(ls)
}
