#' Plotting methods
#'
#' Plotting methods for outputs of \code{bnmonitor} functions
#'
#' @param x an appropriate object
#' @param x The output of node_monitor.
#' @param which select the monitor to plot, either "marginal" or "conditional" (for output of \code{node_monitor} only).
#' @param ... for compatibility
#' @name plot
#' @return A plot specific to the object it is applied to.
NULL


#'@importFrom ggplot2 ggplot xlab ylab theme_minimal ggtitle geom_hline
#'
#' @method plot seq_marg_monitor
#'@export
#' @rdname plot
#'
plot.seq_marg_monitor <- function(x,...){
  temp <- data.frame(x= 1:length(x$Seq_Marg_Monitor), y=x$Seq_Marg_Monitor[1:length(x$Seq_Marg_Monitor)])
  p <- suppressWarnings(ggplot(temp, aes(temp[,1],temp[,2])) + geom_point() +  xlab('Index') + ylab('Standardized Z Statistic') + theme_minimal() + ggtitle(paste0("Marginal Node Monitor for ", x$node.name)) +  geom_hline(yintercept=1.96, linetype="dashed", color = "red") + geom_hline(yintercept=-1.96, linetype="dashed", color = "red"))
  return(p)
}




#' @rdname plot
#'@export
#'
#'@method plot CD

plot.CD <- function(x,...){
  x$plot
}


#'
#'@importFrom ggplot2 ggplot xlab ylab theme_minimal ggtitle geom_hline
#'@rdname plot
#'
#' @method plot seq_cond_monitor
#'@export
#'
#'
plot.seq_cond_monitor <- function(x,...){
  temp <- data.frame(x= 1:length(x$Seq_Cond_Monitor), y=x$Seq_Cond_Monitor[1:length(x$Seq_Cond_Monitor)])
  p <- suppressWarnings(ggplot(temp, aes(temp[,1],temp[,2])) + geom_point() +  xlab('Index') + ylab('Standardized Z Statistic') + theme_minimal() + ggtitle(paste0("Conditional Node Monitor for ", x$node.name)) +  geom_hline(yintercept=1.96, linetype="dashed", color = "red") + geom_hline(yintercept=-1.96, linetype="dashed", color = "red"))
  return(p)
}



#' @importClassesFrom bnlearn bn.fit
#' @importFrom graphics plot.new
#' @importFrom bnlearn arcs
#'@importFrom RColorBrewer brewer.pal
#'@importFrom grDevices colorRampPalette
#'@importFrom qgraph qgraph
#' @method plot node_monitor
#'@export
#'@rdname plot
#'

plot.node_monitor <- function(x, ...){
  nb.cols <- length(names(x$DAG$nodes))
  my.colors <- colorRampPalette(brewer.pal(8, "Blues"))(nb.cols)
  max.val <- ceiling(max(abs(x$Global_Monitor$Score)))
  my.palette <- colorRampPalette(my.colors)(max.val)
  node.colors <- my.palette[floor(abs(x$Global_Monitor$Score))]
  qgraph(x$DAG, color = node.colors, ...)
 # nodes <- create_node_df(n=length(x$DAG$nodes),
 #                         type= names(x$DAG$nodes),
 #                          label=names(x$DAG$nodes),
 #                         style="filled",
 #                          fontcolor="black",
 #                         fillcolor=node.colors, .name_repair = "unique")

 # from.nodes <- arcs(x$DAG)[,1]
 # to.nodes <- arcs(x$DAG)[,2]

 # edges <- create_edge_df(from=match(from.nodes,names(x$DAG$nodes)),
 #                         to=match(to.nodes,names(x$DAG$nodes)))

 # p <- suppressWarnings(create_graph(
 #    nodes_df = nodes,
#    edges_df = edges)
  #  %>%  render_graph( ... )
#  )
  #return(p)
}


#'
#' @method plot influential_obs
#'@export
#'@rdname plot
#'
#' @importFrom ggplot2  xlab ylab theme_minimal

plot.influential_obs <- function(x,...){
  index <- 1:length(x$score)
  value <- x$score
  data <- data.frame(index=index, value = value)
  p <- suppressWarnings(ggplot(data, aes(index, value))+ geom_point() + xlab('Index') + ylab('Leave-One-Out Score') + theme_minimal())
  return(p)
}



#' @rdname plot
#'@export
#'
#'@method plot jeffreys
#'
plot.jeffreys <- function(x,...){
  x$plot
}




#'@export
#' @rdname plot
#'@method plot kl
#'
plot.kl <- function(x,...){
  x$plot
}



#'@importFrom RColorBrewer brewer.pal
#'@importFrom graphics plot.new
#'@importFrom grDevices colorRampPalette
#'@importFrom bnlearn arcs
#'@importFrom qgraph qgraph
#'@method plot final_node_monitor
#'@rdname plot
#'@export
plot.final_node_monitor <- function(x, which, ...){
  if(which!="marginal" & which!="conditional")stop("wrong input for which")
  from.nodes <- arcs(x$DAG)[,1]
  to.nodes <- arcs(x$DAG)[,2]


  #edges <- create_edge_df(from=match(from.nodes,x$Node_Monitor$node),
  #                        to=match(to.nodes,x$Node_Monitor$node))
  l <- length(names(x$DAG$nodes))
  my.colors <-  colorRampPalette(brewer.pal(8, "Greens"))(l)
  max.val <- ceiling(max(abs(x$Node_Monitor$marg.z.score[is.finite(x$Node_Monitor$marg.z.score)])))
  max.val.cond <- ceiling(max(abs(x$Node_Monitor$cond.z.score[is.finite(x$Node_Monitor$cond.z.score)])))
  my.palette <- colorRampPalette(my.colors)(max.val)
  my.palette.cond <- colorRampPalette(my.colors)(max.val.cond)
  node.colors <- my.palette[floor(abs(x$Node_Monitor$marg.z.score))+1]
  node.colors.cond <- my.palette.cond[floor(abs(x$Node_Monitor$cond.z.score))+1]
  #nodes <- create_node_df(n=length(x$Node_Monitor$node),
   #                       type= x$Node_Monitor$node,
  #                        label=x$Node_Monitor$node,
  #                        nodes = x$Node_Monitor$node,
  #                        style="filled",
  #                        fontcolor="black",
  #                        fillcolor=node.colors)

  #nodes.cond <- create_node_df(n=length(x$Node_Monitor$node),
  #                             type= x$Node_Monitor$node,
  #                             label=x$Node_Monitor$node,
  #                             nodes = x$Node_Monitor$node,
  #                             style="filled",
  #                             fontcolor="black",
  #                             fillcolor=node.colors.cond)


  #graph <- create_graph(
  #  nodes_df = nodes,
  #  edges_df = edges)
  #plot <- suppressWarnings(render_graph(graph, title="Marginal Node Monitors", ...))

  #graph.cond <- create_graph(
  #  nodes_df = nodes.cond,
  #  edges_df = edges)
  #plot.cond <- suppressWarnings(render_graph(graph.cond, title="Conditional Node Monitors", ...))
  if(which == "marginal"){qgraph(x$DAG, color = node.colors, ...)
}
  else if(which=="conditional"){ qgraph(x$DAG, color = node.colors.cond, ...)
}
}


#' @importFrom ggplot2 ggtitle xlab ylab theme_minimal theme scale_colour_discrete geom_hline
#' @method plot seq_pa_ch_monitor
#'@export
#'@rdname plot
#' @importFrom ggplot2  xlab ylab theme_minimal

plot.seq_pa_ch_monitor <- function(x,...){
  index <- 1:length(x)
  value <- x[1:length(x)]
  data <- data.frame(index=index, value = value)
  p <- suppressWarnings(ggplot(data, aes(index, value))+ geom_point() + xlab('Relevant sample size') + ylab('Standardized Z Statistic') + theme_minimal() +  geom_hline(yintercept=1.96, linetype="dashed", color = "red") + geom_hline(yintercept=-1.96, linetype="dashed", color = "red"))
  return(p)
}






#'@export
#'@rdname plot
#'@method plot sensitivity
#'
plot.sensitivity <- function(x,...){
  x$plot
}


#'@export
#'@rdname plot
#'@method plot fro
#'
plot.fro <- function(x,...){
  x$plot
}


#' @method plot diameter
#'@export
#'@rdname plot
#' @importFrom qgraph  qgraph
#' @importFrom RColorBrewer brewer.pal
#'@importFrom grDevices colorRampPalette

plot.diameter <- function(x,...){
  nb.cols <- nrow(x$Diameter)
  my.colors <- colorRampPalette(brewer.pal(9, "Reds"))(nb.cols)
  max.val <- 100
  for(i in 1:nrow(x$Diameter)) {if (is.na(x$Diameter[i,2])) {x$Diameter[i,2] <- 0}}
  my.palette <- colorRampPalette(my.colors)(max.val)
  vals <- floor(x$Diameter[,2]*100)
  node.colors <- rep("white",nb.cols)
  for(i in 1:nrow(x$Diameter)) {if(vals[i]!= 0) {node.colors[i] <- my.palette[vals[i]]}}
  qgraph::qgraph(bn.net(x$BN),color=node.colors)
}

#' @method plot edgestrength
#'@export
#'@rdname plot
#' @importFrom qgraph  qgraph

plot.edgestrength <- function(x,...){
  qgraph::qgraph(x[,1:2], edge.labels = round(x[,3],2), edge.width=x[,3]*3,edge.label.cex=1.2)
}

#' @method plot mutualinfo
#'@export
#'@rdname plot
#' @importFrom qgraph  qgraph
#' @importFrom RColorBrewer brewer.pal
#'@importFrom grDevices colorRampPalette
#' @importFrom stats na.omit
#'
plot.mutualinfo <- function(x,...){
  nb.cols <- nrow(x$MutualInfo)
  my.colors <- colorRampPalette(brewer.pal(9, "Reds"))(nb.cols)
  max.val <- (max(na.omit(x$MutualInfo[,2])) + 0.0001)*10^6
  my.palette <- colorRampPalette(my.colors)(max.val)
  vals <- floor(x$MutualInfo[,2]*10^6)
  vals <- ifelse(vals==0,1,vals)
  node.colors <- rep("lightblue",nb.cols)
  for(i in 1:nrow(x$MutualInfo)) {if(!is.na(vals[i])) {node.colors[i] <- my.palette[vals[i]]}}
  qgraph::qgraph(bn.net(x$BN),color=node.colors)
}

#' @method plot dwi
#'@export
#'@rdname plot
#' @importFrom qgraph  qgraph
#' @importFrom RColorBrewer brewer.pal
#'@importFrom grDevices colorRampPalette
#' @importFrom stats na.omit
#'
plot.dwi <- function(x,...){
  nb.cols <- nrow(x$DWI)
  my.colors <- colorRampPalette(brewer.pal(9, "Reds"))(nb.cols)
  max.val <- (max(na.omit(x$DWI[,2])) + 0.0001)*10^6
  my.palette <- colorRampPalette(my.colors)(max.val)
  vals <- floor(x$DWI[,2]*10^6)
  vals <- ifelse(vals==0,1,vals)
  node.colors <- rep("lightblue",nb.cols)
  for(i in 1:nrow(x$DWI)) {if(!is.na(vals[i])) {node.colors[i] <- my.palette[vals[i]]}}
  qgraph::qgraph(x$BN,color=node.colors)
}

#' @method plot ewi
#'@export
#'@rdname plot
#' @importFrom qgraph  qgraph
#' @importFrom RColorBrewer brewer.pal
#'@importFrom grDevices colorRampPalette
#' @importFrom stats na.omit

plot.ewi <- function(x,...){
  nb.cols <- nrow(x$EWI)
  my.colors <- colorRampPalette(brewer.pal(9, "Reds"))(nb.cols)
  max.val <- (max(na.omit(x$EWI[,2])) + 0.0001)*10^6
  my.palette <- colorRampPalette(my.colors)(max.val)
  vals <- floor(x$EWI[,2]*10^6)
  vals <- ifelse(vals==0,1,vals)
  node.colors <- rep("lightblue",nb.cols)
  for(i in 1:nrow(x$EWI)) {if(!is.na(vals[i])) {node.colors[i] <- my.palette[vals[i]]}}
  qgraph::qgraph(bn.net(x$BN),color=node.colors)
}
