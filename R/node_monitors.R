#' Node Monitors
#'
#' Some text here
#'
#' @param dag bnlearn object
#' @param df dataset
#' @param node.name which node
#' @name node
NULL

#'@rdname node
#'@importFrom bnlearn bn.fit as.grain
#'@importFrom gRain querygrain
#'@export
#'
marg.node.monitor <- function(dag,df,node.name){#returns the mth node monitor
  node.idx <- which(colnames(df)==node.name)#TODO test that this exists
  s <- rep(0,dim(df)[1])
  e <- rep(0,dim(df)[1])
  v <- rep(0,dim(df)[1])
  z <- rep(0,dim(df)[1])
  for (i in 1:dim(df)[1]){
    df.cut <- df[1:(i-1),]
    dag.bn.fit <- bn.fit(dag, df.cut)
    dag.grain <- as.grain(dag.bn.fit)
    which.val <- as.numeric(df[i,node.idx])
    querygrain(dag.grain, nodes=colnames(df), type="marginal") ->ev
    p <- unlist(ev[node.name])
    if(any(as.numeric(p)==1) | any(as.numeric(p)==0) ){#| any(as.numeric(p)==1.0) | any(p=="NaN")){
      next
    } else {
      s[i] <- -log(p[which.val])
      e[i] <- -sum((p*log(p)))
      v[i] <- (sum(p*log(p)^2)) -e[i]^2
      z[i] <- (cumsum(s)[i]-cumsum(e)[i])/sqrt(cumsum(v)[i])
    }
  }

  return(z)
}


#'@rdname node
#'@importFrom bnlearn bn.fit as.grain
#'@importFrom gRain querygrain
#'@importFrom purrr map2_dbl
#'@export
marg.node.monitor.graph <- function(dag, df){#node.scores output from global.bn

  num.nodes <- length(nodes(dag))
  dag.bn.fit <- bn.fit(dag, df)
  dag.grain <- as.grain(dag.bn.fit)
  worst.level <- as.numeric(df[dim(df)[1],])
  querygrain(dag.grain, nodes=colnames(df), type="marginal") ->ev
  final.z.score <- ev[match(names(df),names(ev))] %>%  map2_dbl(.y=worst.level,standardize)  #this returns the vvery last marginal

  my.colors = brewer.pal(length(names(dag$nodes)),"Greens")
  max.val <- ceiling(max(abs(final.z.score)))
  my.palette <- colorRampPalette(my.colors)(max.val)
  node.colors <- my.palette[floor(abs(final.z.score))+1]
  nodes <- create_node_df(n=length(names(df)),
                          type= names(df),
                          label=names(df),
                          style="filled",
                          fontcolor="black",
                          fillcolor=node.colors)

  from.nodes <- map(dag$nodes, `[[`, "parents") %>% unlist %>% unname
  to.nodes <- map(dag$nodes, `[[`, "parents") %>% unlist %>% names %>% substr(1,1)

  edges <- create_edge_df(from=match(from.nodes,names(df)),
                          to=match(to.nodes,names(df)))
  create_graph(
    nodes_df = nodes,
    edges_df = edges) %>%
    render_graph(output = "graph",title="Marginal Node Monitors",layout='tree')
}




#'@rdname node
#'@importFrom bnlearn bn.fit as.grain
#'@importFrom gRain setEvidence querygrain
#'@export
cond.node.monitor <- function(dag,df,node.name){#returns the mth node monitor
  node.idx <- which(colnames(df)==node.name)#TODO test that this exists
  s <- rep(0,dim(df)[1]); s.cond <- rep(0,dim(df)[1])
  e <- rep(0,dim(df)[1]); e.cond <- rep(0,dim(df)[1])
  v <- rep(0,dim(df)[1]); v.cond <- rep(0,dim(df)[1])
  z <- rep(0,dim(df)[1]); z.cond <- rep(0,dim(df)[1])
  for (i in 1:dim(df)[1]){
    df.cut <- df[1:(i-1),]
    dag.bn.fit <- bn.fit(dag, df.cut)
    dag.grain <- as.grain(dag.bn.fit)
    dag.ev <- setEvidence(dag.grain,nodes=colnames(df)[-node.idx],states=as.character(unname(unlist(df[i,-node.idx]))))
    p <- unlist(gRain::querygrain(dag.ev, nodes=node.name))
    p.cond <- unlist(querygrain(dag.grain,nodes=node.name))
    which.val <- as.numeric(df[i,node.idx])
    if(any(as.numeric(p)==1) | any(as.numeric(p)==0) | any(as.numeric(p)==1.0) | any(p=="NaN") |
       any(as.numeric(p.cond)==1) | any(as.numeric(p.cond)==0) | any(as.numeric(p.cond)==1.0) | any(p.cond=="NaN")){
      next
    } else {
      s[i] <- -log(p[which.val])
      e[i] <- -sum((p*log(p)))
      v[i] <- (sum(p*log(p)^2)) -e[i]^2
      z[i] <- (cumsum(s)[i]-cumsum(e)[i])/sqrt(cumsum(v)[i])
      s.cond[i] <- -log(p.cond[which.val])
      e.cond[i] <- -sum((p.cond*log(p.cond)))
      v.cond[i] <- (sum(p.cond*log(p.cond)^2)) -e.cond[i]^2
      z.cond[i] <- (cumsum(s.cond)[i]-cumsum(e.cond)[i])/sqrt(cumsum(v.cond)[i])
    }
  }
  z[which(z==0)] <- NA
  z.cond[which(z.cond==0)] <- NA
  return(cbind(z,z.cond))
}

#'@rdname node
#'@importFrom bnlearn bn.fit as.grain
#'@importFrom purrr pmap map2
#'@importFrom RColorBrewer brewer.pal
#'@importFrom grDevices colorRampPalette
#'@importFrom DiagrammeR create_node_df create_edge_df create_graph render_graph
#'@export
cond.node.monitor.graph <- function(dag, df){#node.scores output from global.bn
  dag.bn.fit <- bn.fit(dag, df[1:(dim(df)[1]-1),])
  dag.grain <- as.grain(dag.bn.fit)
  worst.level <- as.numeric(df[dim(df)[1],])
  i_df <- data.frame(
    x=1:length(colnames(df))
  )
  z.scores <- i_df %>% pmap(~pass.ev(.x, df=df, dag.grain=dag.grain)) %>% map2(worst.level,standardize) %>%   unlist %>% unname
  names(z.scores) <- names(df)

  my.colors = brewer.pal(length(names(dag$nodes)),"Greens")
  max.val <- ceiling(max(abs(z.scores)))
  my.palette <- colorRampPalette(my.colors)(max.val)
  node.colors <- my.palette[floor(abs(z.scores))+1]
  nodes <- create_node_df(n=length(names(df)),
                          type= names(df),
                          label=names(df),
                          style="filled",
                          fontcolor="black",
                          fillcolor=node.colors)

  from.nodes <- map(dag$nodes, `[[`, "parents") %>% unlist %>% unname
  to.nodes <-map(dag$nodes, `[[`, "parents") %>% unlist %>% names %>% substr(1,1)

  edges <- create_edge_df(from=match(from.nodes,names(df)),
                          to=match(to.nodes,names(df)))
  create_graph(
    nodes_df = nodes,
    edges_df = edges) %>%
    render_graph(output = "graph",title="Conditional Node Monitors",layout='tree')
}


#'@rdname node
#'@importFrom bnlearn bn.fit as.grain
#'@importFrom gRain querygrain
#'@importFrom tibble new_tibble
#'@export
node.monitor <- function(dag, df){#node.scores output from global.bn
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


  result <-new_tibble(cbind(names(df),as.numeric(marg.z.scores),as.numeric(cond.z.scores)))
  names(result) <- c('node','marg.z.score','cond.z.score')
  return(result)#TODO return the graph as well
}
