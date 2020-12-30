#' Sequential node monitors
#'
#' Sequential marginal and conditional node monitors for a vertex of a Bayesian network.
#'
#' Consider a Bayesian network over variables \eqn{Y_1,\dots,Y_m} and suppose a dataset \eqn{(\boldsymbol{y}_1,\dots,\boldsymbol{y}_n)} has been observed, where \eqn{\boldsymbol{y}_i=(y_{i1},\dots,y_{im})} and \eqn{y_{ij}} is the i-th observation of the j-th variable.
#' Let \eqn{p_i} denote the marginal density of \eqn{Y_j} after the first \eqn{i-1} observations have been processed. Define
#' \deqn{E_i = \sum_{k=1}^Kp_i(d_k)\log(p_i(d_k)),}
#' \deqn{V_i = \sum_{k=1}^K p_i(d_k)\log^2(p_i(d_k))-E_i^2,}
#' where \eqn{(d_1,\dots,d_K)} are the possible values of \eqn{Y_j}. The sequential marginal node monitor for the vertex \eqn{Y_j} is defined as
#' \deqn{Z_{ij}=\frac{-\sum_{k=1}^i\log(p_k(y_{kj}))-\sum_{k=1}^i E_k}{\sqrt{\sum_{k=1}^iV_k}}.}
#' Values of \eqn{Z_{ij}} such that \eqn{|Z_{ij}|> 1.96} can give an indication of a poor model fit for the vertex \eqn{Y_j} after the first i-1 observations have been processed.
#'
#' The sequential conditional node monitor for the vertex \eqn{Y_j} is defined as
#'  \deqn{Z_{ij}=\frac{-\sum_{k=1}^i\log(p_k(y_{kj}|y_{k1},\dots,y_{k(j-1)},y_{k(j+1)},\dots,y_{km}))-\sum_{k=1}^i E_k}{\sqrt{\sum_{k=1}^iV_k}},}
#'  where \eqn{E_k} and \eqn{V_k} are computed with respect to \eqn{p_k(y_{kj}|y_{k1},\dots,y_{k(j-1)},y_{k(j+1)},\dots,y_{km})}. Again, values of \eqn{Z_{ij}} such that \eqn{|Z_{ij}|> 1.96} can give an indication of a poor model fit for the vertex \eqn{Y_j}.
#'
#' @return A vector including the scores \eqn{Z_{ij}}, either marginal or conditional.
#'
#' @param dag an object of class \code{bn} from the \code{bnlearn} package
#' @param df a base R style dataframe
#' @param node.name node over which to compute the monitor
#'
#' @references Cowell, R. G., Dawid, P., Lauritzen, S. L., & Spiegelhalter, D. J. (2006). Probabilistic networks and expert systems: Exact computational methods for Bayesian networks. Springer Science & Business Media.
#' @references Cowell, R. G., Verrall, R. J., & Yoon, Y. K. (2007). Modeling operational risk with Bayesian networks. Journal of Risk and Insurance, 74(4), 795-827.
#'
#' @examples seq_marg_monitor(chds_bn, chds[1:100,], "Events")
#' @examples seq_marg_monitor(chds_bn, chds[1:100,], "Admission")
#'
#' @seealso \code{\link{influential_obs}}, \code{\link{node_monitor}}, \code{\link{seq_node_monitor}}, \code{\link{seq_pa_ch_monitor}}
#' @name seq_node_monitor
NULL

#'@rdname seq_node_monitor
#'@importFrom bnlearn bn.fit as.grain
#'@importFrom gRain querygrain
#'@export
#'
seq_marg_monitor <- function(dag,df,node.name){#returns the mth node monitor
  node.idx <- which(colnames(df)==node.name)#TODO test that this exists
  s <- rep(0,dim(df)[1])
  e <- rep(0,dim(df)[1])
  v <- rep(0,dim(df)[1])
  z <- rep(0,dim(df)[1])
  for (i in 1:dim(df)[1]){
    df.cut <- df[1:i,]
    dag.bn.fit <- suppressWarnings(bn.fit(dag, df.cut))
    dag.grain <- suppressWarnings(as.grain(dag.bn.fit))
    which.val <- suppressWarnings(as.numeric(df[i,node.idx]))
    suppressWarnings(querygrain(dag.grain, nodes=colnames(df), type="marginal")) ->ev
    p <- unlist(ev[node.name])
    if(any(as.numeric(p)==1) | any(as.numeric(p)==0) ){#| any(as.numeric(p)==1.0) | any(p=="NaN")){
      next
    } else {
      s[i] <- suppressWarnings(-log(p[which.val]))
      e[i] <- suppressWarnings(-sum((p*log(p))))
      v[i] <- suppressWarnings((sum(p*log(p)^2)) -e[i]^2)
      z[i] <- (cumsum(s)[i]-cumsum(e)[i])/sqrt(cumsum(v)[i])
    }
  }
  z[which(z==0)] <- NA
  score <- list(Seq_Marg_Monitor = z, node.name = node.name)
  attr(score, 'class') <- 'seq_marg_monitor'
  return(score)
}






#'@rdname seq_node_monitor
#'@importFrom bnlearn bn.fit as.grain
#'@importFrom gRain setEvidence querygrain
#'@export
seq_cond_monitor <- function(dag,df,node.name){#returns the mth node monitor
  node.idx <- which(colnames(df)==node.name)#TODO test that this exists
  s <- rep(0,dim(df)[1]);
  e <- rep(0,dim(df)[1]);
  v <- rep(0,dim(df)[1]);
  z <- rep(0,dim(df)[1]);
  for (i in 1:dim(df)[1]){
    df.cut <- df[1:i,]
    dag.bn.fit <- suppressWarnings(bn.fit(dag, df.cut))
    dag.grain <- suppressWarnings(as.grain(dag.bn.fit))
    dag.ev <- suppressWarnings(setEvidence(dag.grain,nodes=colnames(df)[-node.idx],states=as.character(unname(unlist(df[i,-node.idx])))))
    p <- suppressWarnings(unlist(gRain::querygrain(dag.ev, nodes=node.name)))
    which.val <- as.numeric(df[i,node.idx])
    if(any(as.numeric(p)==1) | any(as.numeric(p)==0) | any(as.numeric(p)==1.0) | any(p=="NaN")){
      next
    } else {
      s[i] <- suppressWarnings(-log(p[which.val]))
      e[i] <- suppressWarnings(-sum((p*log(p))))
      v[i] <- suppressWarnings((sum(p*log(p)^2)) -e[i]^2)
      z[i] <- (cumsum(s)[i]-cumsum(e)[i])/sqrt(cumsum(v)[i])
      }
  }
  z[which(z==0)] <- NA
  score <- list(Seq_Cond_Monitor = z, node.name = node.name)
  attr(score, 'class') <- 'seq_cond_monitor'
  return(score)
}






