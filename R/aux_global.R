
get.pa.combo.score <- function(k, counts.vec, alpha.vec){
  c.vec <- counts.vec[(((k-1)*length(alpha.vec))+1): (((k-1)*length(alpha.vec))+length(alpha.vec))]
  lgamma(sum(alpha.vec)) - lgamma(sum(alpha.vec + c.vec)) + sum(lgamma(alpha.vec + c.vec)) - sum(lgamma(alpha.vec))
}


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
  score <- unlist(-sum(scores.vec))
  return(score)#returns global and pach monitor
}


