#' Node Monitor
#'
#' Some text here
#'
#' @param dag bnlearn object
#' @param df dataset
#' @param plot Boolean
#'
#'@importFrom bnlearn bn.fit as.grain
#'@importFrom gRain querygrain
#'@importFrom gridExtra grid.arrange
#'@importFrom purrr pmap map2 map2_dbl
#'@importFrom RColorBrewer brewer.pal
#'@importFrom grDevices colorRampPalette
#'@importFrom DiagrammeR create_node_df create_edge_df create_graph render_graph
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
 cond.z.scores <- i_df %>% pmap(~pass.ev(.x, df=df, dag.grain=dag.grain)) %>%
    map2(worst.level,standardize) %>%
    unlist %>% unname

  dag.bn.fit.marg <- bn.fit(dag, df)
  dag.grain.marg <- as.grain(dag.bn.fit.marg)
  querygrain(dag.grain.marg, nodes=colnames(df), type="marginal") ->ev
  ev[match(names(df),names(ev))] %>% map2_dbl(.y=worst.level,standardize) -> marg.z.score #this returns the vvery last marginal
  marg.z.scores <- marg.z.score[match(names(df),names(marg.z.score))]


  result <- data.frame(cbind(names(df),as.numeric(as.character(marg.z.scores)),as.numeric(as.character(cond.z.scores))))
  names(result) <- c('node','marg.z.score','cond.z.score')

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
