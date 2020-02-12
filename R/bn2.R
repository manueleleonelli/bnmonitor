#' Integration with \code{bn.fit} objects from \code{bnlearn}
#'
#' Functions that transform an object of class \code{bn.fit} and \code{bn.fit.gnet} (a Gaussian Bayesian network) to objects of class \code{GBN} or \code{CI}.
#'
#'
#'
#'@param bnfit object of class \code{bn.fit}.
#'@return The function \code{bn2gbn} returns an object of class \code{GBN} consisting of a list with entries:
#'\itemize{
#'\item \code{order}: An ordering of the nodes according to the graph.
#'\item \code{mean}: The mean vector of the Gaussian distribution.
#'\item \code{covariance}: The covariance matrix of the Gaussian distribution.
#'}
#'
#'The function \code{bn2ci} returns an object of class \code{CI} consisting of the same list as \code{GBN}, but with the additional entry \code{cond_ind}. \code{cond_ind} is a list where each entry consists of \code{A}, \code{B} and \code{C} corresponding to the conditional independence statements \code{A} independent of \code{B} given \code{C} embedded by the network.
#'@name bn2
NULL

#' @rdname bn2
#' @importFrom bnlearn nodes
#' @importFrom bnlearn node.ordering
#' @importFrom stats coefficients
#'@importClassesFrom bnlearn bn.fit
#'@export

bn2gbn <- function(bnfit){
  B <- matrix(0,length(nodes(bnfit)),length(nodes(bnfit)))
  Phi <- matrix(0,length(nodes(bnfit)),length(nodes(bnfit)))
  b <- rep(0,length(nodes(bnfit)))
  for(i in 1:length(nodes(bnfit))){
    curr <- node.ordering(bnfit)[i]
    b[i] <- coefficients(bnfit)[[curr]][1]
    Phi[i,i] <- bnfit[[curr]]$sd^2
    curr_par <- bnfit[[curr]]$parents
    for (k in 1:length(curr_par)){
      j <- which(node.ordering(bnfit) == curr_par[k])
      B[j,i] <- coefficients(bnfit)[[curr]][k+1]
    }
  }
  cov <- t(solve(diag(rep(1,length(nodes(bnfit))))-B))%*%Phi%*%solve(diag(rep(1,length(nodes(bnfit))))-B)
  mean <- t(solve(diag(rep(1,length(nodes(bnfit))))-B))%*%b
  gbn <- list(order = node.ordering(bnfit), mean = mean, covariance = cov)
  attr(gbn, 'class') <- 'GBN'
  return(gbn)
}

#' @rdname bn2
#' @importFrom bnlearn nodes
#' @importFrom bnlearn node.ordering
#' @importFrom bnlearn arcs
#'@importClassesFrom bnlearn bn.fit
#'@export
#'
bn2ci <- function(bnfit){
  edge <- arcs(bnfit)
  vertex <- node.ordering(bnfit)
  temp <- vector(mode = "list", length = length(vertex) -1)
  edge.temp <- edge
  for (i in 2:length(vertex)){
    A <- vertex[i]
    C <- unique(matrix(edge.temp[ifelse(edge.temp[,2]== vertex[i],T,F),],ncol=2)[,1])
    B <- setdiff(vertex[1:(i-1)],C)
    edge.temp <- matrix(edge.temp[ifelse(edge.temp[,2]== vertex[i],F,T),],ncol=2)
    temp[[i-1]] <- list(A = A, B = B, C = C)
  }
  B <- matrix(0,length(nodes(bnfit)),length(nodes(bnfit)))
  Phi <- matrix(0,length(nodes(bnfit)),length(nodes(bnfit)))
  b <- rep(0,length(nodes(bnfit)))
  for(i in 1:length(nodes(bnfit))){
    curr <- node.ordering(bnfit)[i]
    b[i] <- coefficients(bnfit)[[curr]][1]
    Phi[i,i] <- bnfit[[curr]]$sd^2
    curr_par <- bnfit[[curr]]$parents
    for (k in 1:length(curr_par)){
      j <- which(node.ordering(bnfit) == curr_par[k])
      B[j,i] <- coefficients(bnfit)[[curr]][k+1]
    }
  }
  cov <- t(solve(diag(rep(1,length(nodes(bnfit))))-B))%*%Phi%*%solve(diag(rep(1,length(nodes(bnfit))))-B)
  mean <- t(solve(diag(rep(1,length(nodes(bnfit))))-B))%*%b
  gbn <- list(order = node.ordering(bnfit), mean = mean, covariance = cov, cond_ind = temp)
  attr(gbn,'class') <- 'CI'
  return(gbn)
}
