#' Node monitors
#'
#' Marginal and conditional node monitors for all vertices of a Bayesian network using the full dataset
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
#' @examples node_monitor(chds_bn, chds)
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
node_monitor <- function(dag, df){
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
 attr(result,'class') <- 'node_monitor'
  return(result)
 }

#' Plot for node monitors
#'
#'@importFrom RColorBrewer brewer.pal
#'@importFrom graphics plot.new
#'@importFrom grDevices colorRampPalette
#'@importFrom DiagrammeR create_node_df create_edge_df create_graph render_graph
#'@importFrom bnlearn arcs
#'@param x The output of node_monitor.
#'@param which select the monitor to plot, either "marginal" or "conditional".
#'@param ... additional inputs
#'@method plot node_monitor
#'@export
plot.node_monitor <- function(x, which, ...){
  if(which!="marginal" & which!="conditional")stop("wrong input for which")
  from.nodes <- arcs(x$DAG)[,1]
  to.nodes <- arcs(x$DAG)[,2]

  edges <- create_edge_df(from=match(from.nodes,x$Node_Monitor$node),
                          to=match(to.nodes,x$Node_Monitor$node))

  my.colors = brewer.pal(length(names(x$DAG$nodes)),"Greens")
  max.val <- ceiling(max(abs(x$Node_Monitor$marg.z.score)))
  max.val.cond <- ceiling(max(abs(x$Node_Monitor$cond.z.score)))
  my.palette <- colorRampPalette(my.colors)(max.val)
  my.palette.cond <- colorRampPalette(my.colors)(max.val.cond)
  node.colors <- my.palette[floor(abs(x$Node_Monitor$marg.z.score))+1]
  node.colors.cond <- my.palette.cond[floor(abs(x$Node_Monitor$cond.z.score))+1]
  nodes <- create_node_df(n=length(x$Node_Monitor$node),
                          type= x$Node_Monitor$node,
                          label=x$Node_Monitor$node,
                          nodes = x$Node_Monitor$node,
                          style="filled",
                          fontcolor="black",
                          fillcolor=node.colors)

  nodes.cond <- create_node_df(n=length(x$Node_Monitor$node),
                               type= x$Node_Monitor$node,
                               label=x$Node_Monitor$node,
                               nodes = x$Node_Monitor$node,
                               style="filled",
                               fontcolor="black",
                               fillcolor=node.colors.cond)


  graph <- create_graph(
    nodes_df = nodes,
    edges_df = edges)
  plot <- suppressWarnings(render_graph(graph, title="Marginal Node Monitors",layout = "tree"))

  graph.cond <- create_graph(
    nodes_df = nodes.cond,
    edges_df = edges)
  plot.cond <- suppressWarnings(render_graph(graph.cond, title="Conditional Node Monitors",layout = "tree"))
  if(which == "marginal"){return(plot)}
  else if(which=="conditional"){return(plot.cond)}
}

#' Print of node monitor
#' @importClassesFrom bnlearn bn.fit
#'@export
#'
#'@param x The output of node_monitor
#'@param ... additional inputs
#'
print.node_monitor <- function(x,...){
  print(x$Node_Monitor)
  invisible(x)
}
