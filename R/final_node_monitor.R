#' Final node monitors
#'
#' Marginal and conditional node monitors over the last observation of the data for all vertices of a Bayesian network using the full dataset
#'
#' Consider a Bayesian network over variables \eqn{Y_1,\dots,Y_m} and suppose a dataset \eqn{(\boldsymbol{y}_1,\dots,\boldsymbol{y}_n)} has been observed, where \eqn{\boldsymbol{y}_i=(y_{i1},\dots,y_{im})} and \eqn{y_{ij}} is the i-th observation of the j-th variable.
#' Let \eqn{p_n} denote the marginal density of \eqn{Y_j} after the first \eqn{n-1} observations have been processed. Define
#' \deqn{E_n = \sum_{k=1}^Kp_n(d_k)\log(p_n(d_k)),}
#' \deqn{V_n = \sum_{k=1}^K p_n(d_k)\log^2(p_n(d_k))-E_n^2,}
#' where \eqn{(d_1,\dots,d_K)} are the possible values of \eqn{Y_j}. The marginal node monitor for the vertex \eqn{Y_j} is defined as
#' \deqn{Z_j=\frac{-\log(p_n(y_{ij}))- E_n}{\sqrt{V_n}}.}
#' Higher values of \eqn{Z_j} can give an indication of a poor model fit for the vertex \eqn{Y_j}.
#'
#' The conditional node monitor for the vertex \eqn{Y_j} is defined as
#'  \deqn{Z_j=\frac{-\log(p_n(y_{nj}|y_{n1},\dots,y_{n(j-1)},y_{n(j+1)},\dots,y_{nm}))- E_n}{\sqrt{V_n}},}
#'  where \eqn{E_n} and \eqn{V_n} are computed with respect to \eqn{p_n(y_{nj}|y_{n1},\dots,y_{n(j-1)},y_{n(j+1)},\dots,y_{nm})}. Again, higher values of \eqn{Z_j} can give an indication of a poor model fit for the vertex \eqn{Y_j}.
#'
#' @return A dataframe including the names of the vertices, the marginal node monitors and the conditional node monitors. It also return two plots where vertices with a darker color have a higher marginal z-score or conditional z-score, respectively, in absolute value.
#'
#' @examples final_node_monitor(chds_bn, chds[1:100,])
#' @param dag an object of class \code{bn} from the \code{bnlearn} package
#' @param df a base R style dataframe
#'
#'@importFrom bnlearn bn.fit as.grain
#'@importFrom gRain querygrain
#'@importFrom purrr pmap map2 map2_dbl
#'
#' @references Cowell, R. G., Dawid, P., Lauritzen, S. L., & Spiegelhalter, D. J. (2006). Probabilistic networks and expert systems: Exact computational methods for Bayesian networks. Springer Science & Business Media.
#' @references Cowell, R. G., Verrall, R. J., & Yoon, Y. K. (2007). Modeling operational risk with Bayesian networks. Journal of Risk and Insurance, 74(4), 795-827.
#'
#'
#' @seealso \code{\link{influential_obs}}, \code{\link{node_monitor}}, \code{\link{seq_node_monitor}}, \code{\link{seq_pa_ch_monitor}}
#'@export
final_node_monitor <- function(dag, df){
  #node.scores output from global.bn
  num.nodes <- length(nodes(dag))
  dag.bn.fit <- bn.fit(dag, df[1:(dim(df)[1]-1),])
  dag.grain <- as.grain(dag.bn.fit)
  worst.level <- as.numeric(df[dim(df)[1],])
  i_df <- data.frame(
    x=1:length(colnames(df))
  )
 cond.z.scores <- as.numeric(as.character(i_df %>% pmap(~pass.ev(.x, df=df, dag.grain=dag.grain)) %>% map2(worst.level,standardize) %>% unlist %>% unname))

  dag.bn.fit.marg <- bn.fit(dag, df)
  dag.grain.marg <- as.grain(dag.bn.fit.marg)
  querygrain(dag.grain.marg, nodes=colnames(df), type="marginal") ->ev
  ev[match(names(df),names(ev))] %>% map2_dbl(.y=worst.level,standardize) -> marg.z.score #this returns the vvery last marginal
  marg.z.scores <- as.numeric(as.character(marg.z.score[match(names(df),names(marg.z.score))]))



 temp <- data.frame(node = names(df),marg.z.score = marg.z.scores, cond.z.score = cond.z.scores)
 result <- list(Node_Monitor = temp, DAG = dag)
 attr(result,'class') <- 'final_node_monitor'
  return(result)
 }

