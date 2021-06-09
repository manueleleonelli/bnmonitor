
#' Sequential parent-child node monitors
#'
#' Sequential node monitor for a vertex of a Bayesian network for a specific configuration of its parents
#'
#' Consider a Bayesian network over variables \eqn{Y_1,\dots,Y_m} and suppose a dataset \eqn{(\boldsymbol{y}_1,\dots,\boldsymbol{y}_n)} has been observed, where \eqn{\boldsymbol{y}_i=(y_{i1},\dots,y_{im})} and \eqn{y_{ij}} is the i-th observation of the j-th variable.
#' Consider a configuration \eqn{\pi_j} of the parents and consider the sub-vector \eqn{\boldsymbol{y}'=(\boldsymbol{y}_1',\dots,\boldsymbol{y}_{N'}')} of \eqn{(\boldsymbol{y}_1,\dots,\boldsymbol{y}_n)} including observations where the parents of \eqn{Y_j} take value \eqn{\pi_j} only.
#' Let \eqn{p_i(\cdot|\pi_j)} be the conditional distribution of \eqn{Y_j} given that its parents take value \eqn{\pi_j} after the first i-1 observations have been processed. Define
#' \deqn{E_i = \sum_{k=1}^Kp_i(d_k|\pi_j)\log(p_i(d_k|\pi_j)),}
#' \deqn{V_i = \sum_{k=1}^K p_i(d_k|\pi_j)\log^2(p_i(d_k|\pi_j))-E_i^2,}
#' where \eqn{(d_1,\dots,d_K)} are the possible values of \eqn{Y_j}. The sequential parent-child node monitor for the vertex \eqn{Y_j} and parent configuration \eqn{\pi_j} is defined as
#'  \deqn{Z_{ij}=\frac{-\sum_{k=1}^i\log(p_k(y_{kj}'|\pi_j))-\sum_{k=1}^i E_k}{\sqrt{\sum_{k=1}^iV_k}}.}
#'  Values of \eqn{Z_{ij}} such that \eqn{|Z_{ij}|> 1.96} can give an indication of a poor model fit for the vertex \eqn{Y_j} after the first i-1 observations have been processed.
#'
#' @return A vector including the scores \eqn{Z_{ij}}.
#'
#'@param dag an object of class \code{bn} from the \code{bnlearn} package
#'@param df a base R style dataframe
#'@param node.name node over which to compute the monitor
#'@param pa.names vector including the names of the parents of \code{node.name}
#'@param pa.val vector including the levels of \code{pa.names} considered
#'@param alpha single integer. By default, the number of max levels in \code{df}
#'
#'
#' @examples seq_pa_ch_monitor(chds_bn, chds, "Events", "Social", "High", 3)
#'
#'@importFrom purrr map_dbl map map_int
#'@importFrom rlang sym
#' @importFrom dplyr "%>%"
#'@importFrom dplyr if_else filter
#' @importFrom reshape2 melt
#'
#' @references Cowell, R. G., Dawid, P., Lauritzen, S. L., & Spiegelhalter, D. J. (2006). Probabilistic networks and expert systems: Exact computational methods for Bayesian networks. Springer Science & Business Media.
#' @references Cowell, R. G., Verrall, R. J., & Yoon, Y. K. (2007). Modeling operational risk with Bayesian networks. Journal of Risk and Insurance, 74(4), 795-827.
#'
#'
#'
#' @seealso \code{\link{influential_obs}}, \code{\link{node_monitor}}, \code{\link{seq_node_monitor}}, \code{\link{seq_pa_ch_monitor}}
#'@export

seq_pa_ch_monitor <- function(dag, df, node.name, pa.names, pa.val, alpha = "default"){#takes input from bnlearn
  if (alpha == "default") alpha <- max(sapply(df, nlevels))
  for (i in length(pa.names)){
    df <- dplyr::filter(df,!!(sym(pa.names[i]))==pa.val[i])
  }
  nodes <-nodes(dag)
  num.nodes <- length(dag$nodes)
  num.ch <- map(dag$nodes, `[[`, "children")  %>% map_int(length)
  num.values <- map_int(1:num.nodes, function(i){length(unique(df[,i]))})
  node.idx <- which(colnames(df)==node.name)
  alpha.vec <- rep(alpha/num.values[node.idx], num.values[node.idx])
  #alpha.vec <- c(0.8,0.2)
  new.alpha <- alpha.vec #this is the one with learning
  counts <- rep(0,num.values[node.idx])
  #idx <- which(all(df[,which(colnames(df) %in% pa.names)]==pa.val))

  #c <- num.values[node.idx]
  s <- rep(0,dim(df)[1]) ;s.learn <- rep(0,dim(df)[1])
  e <- rep(0,dim(df)[1]) ;e.learn <- rep(0,dim(df)[1])
  v <- rep(0,dim(df)[1]) ;v.learn <- rep(0,dim(df)[1])
  z <- rep(0,dim(df)[1]) ;z.learn <- rep(0,dim(df)[1])
  n <- 0

  for(i in 1:(dim(df)[1]-1)){
    if(all(df[i,pa.names]==pa.val)){
      counts[as.numeric(df[i,node.idx])] <- counts[as.numeric(df[i,node.idx])] + 1
      new.alpha[as.numeric(df[i,node.idx])] <- new.alpha[as.numeric(df[i,node.idx])] + 1
      which.val <- as.numeric(df[(i+1),node.idx])
      p <- (alpha.vec)/sum(alpha.vec)
      n <- n+1

      s[i] <- -log(p)[which.val]
      e[i] <- -sum((p*log(p)))

      p.learn <- (new.alpha) /sum(new.alpha)
      #p.learn <- lgamma(sum(new.alpha))-sum(lgamma(new.alpha)) + sum(lgamma(new.alpha + counts)) - lgamma(sum(new.alpha + counts))
      s.learn[i] <- -log(p.learn)[which.val]

      e.learn[i] <- -sum((p.learn*log(p.learn)))
      pi <- e[i]/n; pi.learn <- e.learn[i]/n;

      v[i] <- (sum(p*log(p)^2)) -e[i]^2#+(pi^2/12)
      #v[i] <- n*pi*(1-pi)
      z[i] <- (cumsum(s)[i]-cumsum(e)[i])/sqrt(cumsum(v)[i])
      v.learn[i] <- (sum(p.learn*log(p.learn)^2)) -e.learn[i]^2#-(c^2/12)
      z.learn[i] <- (cumsum(s.learn)[i]-cumsum(e.learn)[i])/sqrt(cumsum(v.learn)[i])
    } else {
      new.alpha[as.numeric(df[i,node.idx])] <- new.alpha[as.numeric(df[i,node.idx])] + 1
      next
    }
  }

  z.learn[which(z.learn==0)] <- NA
  score <-  z.learn
  attr(score, 'class') <- 'seq_pa_ch_monitor'
  return(score)
}




