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
#' @param df a base R style dataframe
#' @param alpha single integer, usually the number of max levels in \code{df}
#' @param plot boolean value. If \code{TRUE} the function returns a plot.
#'
#' @importClassesFrom bnlearn bn.fit
#'@importFrom purrr  map_dbl
#'@importFrom ggplot2 ggplot xlab ylab theme_minimal ggtitle
#'
#' @examples influential_obs(chds_bn, chds[1:100,], 3, FALSE)
#'
#' @seealso \code{\link{influential_obs}}, \code{\link{node_monitor}}, \code{\link{seq_node_monitor}}, \code{\link{seq_pa_ch_monitor}}
#'@export
#'
influential_obs <- function(dag, df, alpha, plot = TRUE){#j is the index of the parent set
  result <- rep(0,nrow(df))
  for (i in 1:nrow(df)){
    result[i] <- sum(as.numeric(as.character(map_dbl(.x=1:length(dag$nodes), dag, alpha, df[-i,], .f= global.monitor.bn.node))))
  }
  total <- sum(as.numeric(as.character(map_dbl(.x=1:length(dag$nodes), dag, alpha, df, .f= global.monitor.bn.node))))
  score <- abs(total - result)
 if(plot == TRUE){
   score <- data.frame(score = score)
    p <- suppressWarnings(ggplot(score, aes(x = 1:nrow(df), y = score)) + geom_point() +  xlab('Index') + ylab('Leave-One-Out Score') + theme_minimal() )
  return(list(score = score, plot = p ))
 }
  return(score)#returns global and pach monitor
}