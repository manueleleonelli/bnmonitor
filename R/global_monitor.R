#' Global monitor
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
#' @param alpha single integer, usually the number of max levels in \code{df}
#' @param plot boolean value. If \code{TRUE} the function returns a plot
#'
#' @examples global_monitor(chds_bn, chds, 3, FALSE)
#'
#' @importClassesFrom bnlearn bn.fit
#'@importFrom purrr map map_int map_dbl
#'@importFrom rlang is_empty syms
#'@importFrom dplyr "%>%"
#'@importFrom dplyr count
#'@importFrom tidyr complete
#'@importFrom dplyr if_else
#'@importFrom dplyr pull
#'@importFrom RColorBrewer brewer.pal
#'@importFrom grDevices colorRampPalette
#'@importFrom DiagrammeR create_node_df create_edge_df create_graph render_graph
#'
#'
#' @references Cowell, R. G., Dawid, P., Lauritzen, S. L., & Spiegelhalter, D. J. (2006). Probabilistic networks and expert systems: Exact computational methods for Bayesian networks. Springer Science & Business Media.
#' @references Cowell, R. G., Verrall, R. J., & Yoon, Y. K. (2007). Modeling operational risk with Bayesian networks. Journal of Risk and Insurance, 74(4), 795-827.
#'
#' @seealso \code{\link{influential_obs}}, \code{\link{node_monitor}}, \code{\link{seq_node_monitor}}, \code{\link{seq_pa_ch_monitor}}
#'
#'@export
#'
global_monitor <- function(dag, df, alpha, plot = TRUE){
  node.scores <- as.numeric(as.character(map_dbl(.x=1:length(dag$nodes), dag, alpha, df, .f= global.monitor.bn.node)))
  result <- data.frame(Vertex = names(dag$nodes), Score = node.scores)
  result <- list(Global_Monitor = result, DAG = dag)
  attr(result, 'class') <- 'global'
  return(result)
}

# Plot of global monitor

plot.global <- function(result){
    my.colors = brewer.pal(length(names(result$DAG$nodes)),"Blues")
    max.val <- ceiling(max(abs(result$Global_Monitor$Score)))
    my.palette <- colorRampPalette(my.colors)(max.val)
    node.colors <- my.palette[floor(abs(result$Global_Monitor$Score))]
    nodes <- create_node_df(n=length(result$DAG$nodes),
                            type= names(result$DAG$nodes),
                            label=names(result$DAG$nodes),
                            style="filled",
                            fontcolor="black",
                            fillcolor=node.colors, .name_repair = "unique")

    from.nodes <- arcs(result$DAG)[,1]
    to.nodes <- arcs(result$DAG)[,2]

    edges <- create_edge_df(from=match(from.nodes,names(result$DAG)),
                            to=match(to.nodes,names(result$DAG)))

    p <- suppressWarnings(create_graph(
      nodes_df = nodes,
      edges_df = edges) %>%
        render_graph(title="Global Monitors",layout="tree"))
    return(p)
}
