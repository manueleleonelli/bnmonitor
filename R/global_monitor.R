#' Global monitor
#'
#' Negative marginal log-likelihood of the model
#'
#'
#' @return A numerical value
#'
#' @param dag an object of class \code{bn} from the \code{bnlearn} package
#' @param df a base R style dataframe
#' @param alpha single integer. By default, number of max levels in \code{df}
#'
#' @examples global_monitor(chds_bn, chds, 3)
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
#' @seealso \code{\link{node_monitor}}, \code{\link{influential_obs}}, \code{\link{final_node_monitor}}, \code{\link{seq_node_monitor}}, \code{\link{seq_pa_ch_monitor}}
#'
#'@export
#'
global_monitor <- function(dag, df, alpha = "default"){
  if (alpha == "default") alpha <- max(sapply(df, nlevels))
  node.scores <- as.numeric(as.character(map_dbl(.x=1:length(dag$nodes), dag, alpha, df, .f= global.monitor.bn.node)))
  result <- data.frame(Vertex = names(dag$nodes), Score = node.scores)
#  result <- list(Global_Monitor = result, DAG = dag)
#  attr(result, 'class') <- 'node_monitor'
  return(sum(result$Score))
}


