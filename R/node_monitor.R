#' Node monitors
#'
#' Marginal and conditional node monitors for all vertices of a Bayesian network using the full dataset
#'
#' Consider a Bayesian network over variables \eqn{Y_1,\dots,Y_m} and suppose a dataset \eqn{(\boldsymbol{y}_1,\dots,\boldsymbol{y}_n)} has been observed, where \eqn{\boldsymbol{y}_i=(y_{i1},\dots,y_{im})} and \eqn{y_{ij}} is the i-th observation of the j-th variable.
#' Let \eqn{p_i} denote the marginal density of \eqn{Y_j} after the first \eqn{i-1} observations have been processed. Define
#' \deqn{E_i = \sum_{k=1}^Kp_i(d_k)\log(p_i(d_k)),}
#' \deqn{V_i = \sum_{k=1}^K p_i(d_k)\log^2(p_i(d_k))-E_i^2,}
#' where \eqn{(d_1,\dots,d_K)} are the possible values of \eqn{Y_j}. The marginal node monitor for the vertex \eqn{Y_j} is defined as
#' \deqn{Z_j=\frac{-\sum_{i=1}^N\log(p_i(y_{ij}))-\sum_{i=1}^N E_i}{\sqrt{\sum_{i=1}^NV_i}}.}
#' Values of \eqn{Z_j} such that \eqn{|Z_j|> 1.96} can give an indication of a poor model fit for the vertex \eqn{Y_j}.
#'
#' The conditional node monitor for the vertex \eqn{Y_j} is defined as
#'  \deqn{Z_j=\frac{-\sum_{i=1}^N\log(p_i(y_{ij}|y_{i1},\dots,y_{i(j-1)},y_{i(j+1)},\dots,y_{im}))-\sum_{i=1}^N E_i}{\sqrt{\sum_{i=1}^NV_i}},}
#'  where \eqn{E_i} and \eqn{V_i} are computed with respect to \eqn{p_i(y_{ij}|y_{i1},\dots,y_{i(j-1)},y_{i(j+1)},\dots,y_{im})}. Again, values of \eqn{Z_j} such that \eqn{|Z_j|> 1.96} can give an indication of a poor model fit for the vertex \eqn{Y_j}.
#'
#' @return A dataframe including the names of the vertices, the marginal node monitors and the conditional node monitors. It also return two plots where vertices with a darker color have a higher marginal z-score or conditional z-score, respectively, in absolute value.
#'
#' @param dag an object of class \code{bn} from the \code{bnlearn} package
#' @param df a base R style dataframe
#' @param plot boolean value. If \code{TRUE} the function returns a plot.
#'
#'@importFrom bnlearn bn.fit as.grain
#'@importFrom gRain querygrain
#'@importFrom gridExtra grid.arrange
#'@importFrom purrr pmap map2 map2_dbl
#'@importFrom RColorBrewer brewer.pal
#'@importFrom grDevices colorRampPalette
#'@importFrom DiagrammeR create_node_df create_edge_df create_graph render_graph
#'
#' @references Cowell, R. G., Dawid, P., Lauritzen, S. L., & Spiegelhalter, D. J. (2006). Probabilistic networks and expert systems: Exact computational methods for Bayesian networks. Springer Science & Business Media.
#' @references Cowell, R. G., Verrall, R. J., & Yoon, Y. K. (2007). Modeling operational risk with Bayesian networks. Journal of Risk and Insurance, 74(4), 795-827.
#'
#'
#' @seealso \code{\link{influential_obs}}, \code{\link{node_monitor}}, \code{\link{seq_node_monitor}}, \code{\link{seq_pa_ch_monitor}}
#'@export
node_monitor <- function(dag, df , plot = TRUE){
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


  result <- data.frame(node = names(df),marg.z.score = marg.z.scores, cond.z.score = cond.z.scores)

  if(plot == TRUE){
    from.nodes <- map(dag$nodes, `[[`, "parents") %>% unlist %>% unname
    to.nodes <- map(dag$nodes, `[[`, "parents") %>% unlist %>% names %>% substr(1,1)

    edges <- create_edge_df(from=match(from.nodes,names(df)),
                            to=match(to.nodes,names(df)))

  my.colors = brewer.pal(length(names(dag$nodes)),"Greens")
  max.val <- ceiling(max(abs(marg.z.scores)))
  max.val.cond <- ceiling(max(abs(cond.z.scores)))
  my.palette <- colorRampPalette(my.colors)(max.val)
  my.palette.cond <- colorRampPalette(my.colors)(max.val.cond)
  node.colors <- my.palette[floor(abs(marg.z.scores))+1]
  node.colors.cond <- my.palette.cond[floor(abs(cond.z.scores))+1]
  nodes <- create_node_df(n=length(names(df)),
                          type= names(df),
                          label=names(df),
                          nodes = names(df),
                          style="filled",
                          fontcolor="black",
                          fillcolor=node.colors)

  nodes.cond <- create_node_df(n=length(names(df)),
                          type= names(df),
                          label=names(df),
                          nodes = names(df),
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
    return(list(table =result, marg.plot = plot, cond.plot = plot.cond))#TODO return the graph as well

    }
 return(result)
 }
