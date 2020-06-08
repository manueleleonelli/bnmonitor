#' Global monitors and the dependent functions

#' @param node.idx index from the dag nodes
#' @param dag bnlearn object specifying a dag
#' @param alpha single integer, usually the number of max levels in df
#' @param df a base R style dataframe
#'
#' @name global
NULL


#'@rdname global
#' @importClassesFrom bnlearn bn.fit
#'@importFrom purrr map map_int map_dbl
#'@importFrom rlang is_empty syms
#'@importFrom dplyr "%>%"
#'@importFrom dplyr count
#'@importFrom tidyr complete
#'@importFrom dplyr if_else
#'@importFrom dplyr pull
#'@export

global.monitor.bn.node <- function(node.idx,dag,alpha,df){#j is the index of the parent set
  num.nodes <- length(dag$nodes)
  pa.val <- map(dag$nodes, `[[`, "parents") #proeprocessing
  num.values <- map_int(1:num.nodes, function(i){length(unique(df[,i]))})
  num.pa.combo <- if_else(is_empty(pa.val[node.idx]),1,prod(num.values[which(colnames(df) %in% unlist(pa.val[node.idx]))]))


  alpha.vec <- rep(alpha/num.values[node.idx], num.values[node.idx])
  pa.names <-c(unlist(pa.val[node.idx] ,use.names = FALSE), names(pa.val)[[node.idx]])
  n <- 0
  counts.vec <- df %>% count(!!!(syms(pa.names))) %>% complete(!!!(syms(pa.names)),fill = list(n = 0)) %>% pull(n)
  scores.vec <- map_dbl(1:num.pa.combo, ~get.pa.combo.score(.x, counts.vec, alpha.vec))
  score <- unlist(sum(scores.vec))
  return(score)#returns global and pach monitor
}

#'@rdname global
#'@importFrom purrr map_dbl
#'@export
#'
global.monitor <- function(dag, alpha, df){#node.scores output from global.bn

  node.scores <- map_dbl(.x=1:length(dag$nodes), dag, alpha, df, .f= global.monitor.bn.node)
  result <- data.frame(Vertex = cbind(names(dag$nodes), Score = node.scores))
  return(result)
}


#'@rdname global
#'@importFrom RColorBrewer brewer.pal
#'@importFrom purrr map_dbl map
#'@importFrom grDevices colorRampPalette
#'@importFrom DiagrammeR create_node_df create_edge_df create_graph render_graph
#'
#'@export
#'
global.monitor.graph <- function(dag, alpha, df){#node.scores output from global.bn

  node.scores <- map_dbl(.x=1:length(dag$nodes), dag, alpha, df, .f= global.monitor.bn.node)
  result <- data.frame(Vertex = cbind(names(dag$nodes),Score = node.scores))


  my.colors = brewer.pal(length(names(dag$nodes)),"Blues")
  max.val <- ceiling(max(abs(node.scores)))
  my.palette <- colorRampPalette(my.colors)(max.val)
  node.colors <- my.palette[floor(abs(node.scores))]
  nodes <- create_node_df(n=length(dag$nodes),
                          type= names(dag$nodes),
                          label=names(dag$nodes),
                          style="filled",
                          fontcolor="black",
                          fillcolor=node.colors, .name_repair = "unique")

  from.nodes <- map(dag$nodes, `[[`, "parents") %>% unlist %>% unname
  to.nodes <-map(dag$nodes, `[[`, "parents") %>% unlist %>% names %>% substr(1,1)

  edges <- create_edge_df(from=match(from.nodes,names(dag$nodes)),
                          to=match(to.nodes,names(dag$nodes)))
  create_graph(
    nodes_df = nodes,
    edges_df = edges) %>%
    render_graph(title="Global Monitors",layout="tree")#remove tree when you've got model 0 to model 1 diff
}
