#' Influential observations
#'
#' Influence of a single observation to the global monitor
#'
#' Consider a Bayesian network over variables \eqn{Y_1,\dots,Y_m} and suppose a dataset \eqn{(\boldsymbol{y}_1,\dots,\boldsymbol{y}_n)} has been observed, where \eqn{\boldsymbol{y}_i=(y_{i1},\dots,y_{im})} and \eqn{y_{ij}} is the i-th observation of the j-th variable. Define \eqn{\boldsymbol{y}_{-i}=(\boldsymbol{y}_1,\dots,\boldsymbol{y}_{i-1},\boldsymbol{y}_{i+1},\dots,\boldsymbol{y}_n)}.
#' The influence of an observation to the global monitor is defined as
#' \deqn{|\log(p(\boldsymbol{y}_1,\dots,\boldsymbol{y}_n)) - \log(p(\boldsymbol{y}_{-i}))|.}
#' High values of this index denote observations that highly contribute to the likelihood of the model.
#'
#' @return A vector including the influence of each observation.
#'
#' @param dag an object of class \code{bn} from the \code{bnlearn} package
#' @param data a base R style dataframe
#' @param alpha single integer. By default, the number of max levels in \code{data}
#'
#' @importClassesFrom bnlearn bn.fit
#'@importFrom purrr  map_dbl
#'@importFrom ggplot2 ggplot xlab ylab theme_minimal ggtitle
#'
#' @examples influential_obs(chds_bn, chds[1:100,], 3)
#'
#' @seealso \code{\link{influential_obs}}, \code{\link{node_monitor}}, \code{\link{seq_node_monitor}}, \code{\link{seq_pa_ch_monitor}}
#'@export
#'
influential_obs <- function(dag,data,alpha = "default"){
  if (alpha == "default") alpha <- max(sapply(data, nlevels))
  un <- unique(data)
  result <- rep(0,nrow(un))
  total <- global_monitor(dag,data)
  for(i in 1:nrow(un)){
    result[i] <- global_monitor(dag,rbind(data,un[i,]),alpha)
  }
  score <- result- total
  out <- data.frame(un,score)
  attr(out,'class') <- c('influential_obs','data.frame')

  return(out)
}
