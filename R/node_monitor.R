#' Node monitor
#'
#' Contribution of each vertex of a Bayesian network to the global monitor
#'
#' Consider a Bayesian network over variables \eqn{Y_1,\dots,Y_m} and suppose a dataset \eqn{(\boldsymbol{y}_1,\dots,\boldsymbol{y}_n)} has been observed, where \eqn{\boldsymbol{y}_i=(y_{i1},\dots,y_{im})} and \eqn{y_{ij}} is the i-th observation of the j-th variable. The global monitor is defined as the negative log-likelihood of the model, i.e.
#' \deqn{-\log(p(\boldsymbol{y}_1,\dots,\boldsymbol{y}_n))= - \sum_{j=1}^m\sum_{i=1}^n \log(p(y_{ij} | \pi_{ij})),}
#' where \eqn{\pi_{ij}} is the value of the parents of \eqn{Y_j} for the i-th observation. The contribution of the j-th vertex to the global monitor is thus
#' \deqn{-\sum_{i=1}^n\log(p(y_{ij}|\pi_{ij})).}
#'
#' @return A dataframe including the name of the vertices and the contribution of the vertices to the global monitor. It also returns a plot where vertices with higher contributions in absolute value are darker.
#'
#' @param dag an object of class \code{bn} from the \code{bnlearn} package
#' @param df a base R style dataframe
#' @param alpha single integer. By default, number of max levels in \code{df}
#'
#' @examples node_monitor(chds_bn, chds, 3)
#'
#' @importClassesFrom bnlearn bn.fit
#'@importFrom purrr map map_int map_dbl
#'@importFrom rlang is_empty syms
#'@importFrom dplyr "%>%"
#'@importFrom dplyr count
#'@importFrom tidyr complete
#'@importFrom dplyr if_else
#'@importFrom dplyr pull
#'@importFrom graphics plot.default
#'
#'
#' @references Cowell, R. G., Dawid, P., Lauritzen, S. L., & Spiegelhalter, D. J. (2006). Probabilistic networks and expert systems: Exact computational methods for Bayesian networks. Springer Science & Business Media.
#' @references Cowell, R. G., Verrall, R. J., & Yoon, Y. K. (2007). Modeling operational risk with Bayesian networks. Journal of Risk and Insurance, 74(4), 795-827.
#'
#' @seealso \code{\link{global_monitor}}, \code{\link{influential_obs}}, \code{\link{final_node_monitor}}, \code{\link{seq_node_monitor}}, \code{\link{seq_pa_ch_monitor}}
#'
#'@export
#'
node_monitor <- function(dag, df, alpha = "default"){
  if (alpha == "default") alpha <- max(sapply(df, nlevels))
  node.scores <- as.numeric(as.character(map_dbl(.x=1:length(dag$nodes), dag, alpha, df, .f= global.monitor.bn.node)))
  result <- data.frame(Vertex = names(dag$nodes), Score = node.scores)
  result <- list(Global_Monitor = result, DAG = dag)
  attr(result, 'class') <- 'node_monitor'
  return(result)
}


