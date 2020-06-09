#' Sequential Node Monitors
#'
#' Some text here
#'
#' @param dag bnlearn object
#' @param df dataset
#' @param node.name which node
#' @param plot Boolean
#' @name seq_node_monitor
NULL

#'@rdname seq_node_monitor
#'@importFrom bnlearn bn.fit as.grain
#'@importFrom gRain querygrain
#'@importFrom ggplot2 ggplot xlab ylab theme_minimal ggtitle
#'@export
#'
seq_marg_monitor <- function(dag,df,node.name, plot = TRUE){#returns the mth node monitor
  node.idx <- which(colnames(df)==node.name)#TODO test that this exists
  s <- rep(0,dim(df)[1])
  e <- rep(0,dim(df)[1])
  v <- rep(0,dim(df)[1])
  z <- rep(0,dim(df)[1])
  for (i in 1:dim(df)[1]){
    df.cut <- df[1:i,]
    dag.bn.fit <- suppressWarnings(bn.fit(dag, df.cut))
    dag.grain <- suppressWarnings(as.grain(dag.bn.fit))
    which.val <- suppressWarnings(as.numeric(df[i,node.idx]))
    suppressWarnings(querygrain(dag.grain, nodes=colnames(df), type="marginal")) ->ev
    p <- unlist(ev[node.name])
    if(any(as.numeric(p)==1) | any(as.numeric(p)==0) ){#| any(as.numeric(p)==1.0) | any(p=="NaN")){
      next
    } else {
      s[i] <- suppressWarnings(-log(p[which.val]))
      e[i] <- suppressWarnings(-sum((p*log(p))))
      v[i] <- suppressWarnings((sum(p*log(p)^2)) -e[i]^2)
      z[i] <- (cumsum(s)[i]-cumsum(e)[i])/sqrt(cumsum(v)[i])
    }
  }
  z[which(z==0)] <- NA
  score <- data.frame(z.score = z)
  if(plot == TRUE){
  z.score <- z
  p <- suppressWarnings(ggplot(score, aes(x = 1:nrow(df), y = z.score)) + geom_point() +  xlab('Index') + ylab('Standardized Z Statistic') + theme_minimal() + ggtitle(paste0("Marginal Node Monitor for ",node.name)))
  return(list(score = score, plot =p ))
  }
  return(score = score)
}




#'@rdname seq_node_monitor
#'@importFrom bnlearn bn.fit as.grain
#'@importFrom gRain setEvidence querygrain
#'@importFrom ggplot2 ggplot xlab ylab theme_minimal ggtitle
#'@export
seq_cond_monitor <- function(dag,df,node.name, plot = TRUE){#returns the mth node monitor
  node.idx <- which(colnames(df)==node.name)#TODO test that this exists
  s <- rep(0,dim(df)[1]);
  e <- rep(0,dim(df)[1]);
  v <- rep(0,dim(df)[1]);
  z <- rep(0,dim(df)[1]);
  for (i in 1:dim(df)[1]){
    df.cut <- df[1:i,]
    dag.bn.fit <- suppressWarnings(bn.fit(dag, df.cut))
    dag.grain <- suppressWarnings(as.grain(dag.bn.fit))
    dag.ev <- suppressWarnings(setEvidence(dag.grain,nodes=colnames(df)[-node.idx],states=as.character(unname(unlist(df[i,-node.idx])))))
    p <- suppressWarnings(unlist(gRain::querygrain(dag.ev, nodes=node.name)))
    which.val <- as.numeric(df[i,node.idx])
    if(any(as.numeric(p)==1) | any(as.numeric(p)==0) | any(as.numeric(p)==1.0) | any(p=="NaN")){
      next
    } else {
      s[i] <- suppressWarnings(-log(p[which.val]))
      e[i] <- suppressWarnings(-sum((p*log(p))))
      v[i] <- suppressWarnings((sum(p*log(p)^2)) -e[i]^2)
      z[i] <- (cumsum(s)[i]-cumsum(e)[i])/sqrt(cumsum(v)[i])
      }
  }
  z[which(z==0)] <- NA
  score <- data.frame(z.score = z)
  if(plot == TRUE){
    z.score <- z
    p <- suppressWarnings(ggplot(score, aes(x = 1:nrow(df), y = z.score)) + geom_point() +  xlab('Index') + ylab('Standardized Z Statistic') + theme_minimal() + ggtitle(paste0("Conditional Node Monitor for ",node.name)))
    return(list(score = score, plot =p ))
  }
  return(score = score)

}
