
#' the sequential parent child node monitor for Bayesian networks
#'@param dframe base R style data frame
#'@param dag bnlearn object with DAG structure, not a fitte dobject
#'@param node.name name of child node
#'@param pa.names names of parents in the network
#'@param pa.val names of values to test
#'@param alpha prior
#' @name  sequential
NULL


#' @rdname sequential
#'@importFrom purrr map_dbl map map_int
#'@importFrom rlang sym
#' @importFrom dplyr "%>%"
#'@importFrom dplyr if_else filter
#'@importFrom graphics plot
#'@export

sqntl.pa.ch.monitor <- function(dframe, dag, node.name, pa.names, pa.val,alpha){#takes input from bnlearn
  df <- dplyr::filter(dframe,!!(sym(pa.names))==pa.val)
  nodes <-nodes(dag)
  num.nodes <- length(dag$nodes)
  num.ch <- map(dag$nodes, `[[`, "children")  %>% map_int(length)
  num.values <- map_int(1:num.nodes, function(i){length(unique(dframe[,i]))})
  node.idx <- which(colnames(df)==node.name)
  alpha.vec <- rep(alpha/num.values[node.idx], num.values[node.idx])
  #alpha.vec <- c(0.8,0.2)
  new.alpha <- alpha.vec #this is the one with learning
  counts <- rep(0,num.values[node.idx])
  #idx <- which(all(df[,which(colnames(df) %in% pa.names)]==pa.val))

  #c <- num.values[node.idx]
  s <- rep(0,dim(df)[1]) ;s.learn <- rep(0,dim(df)[1])
  e <- rep(0,dim(df)[1]) ;e.learn <- rep(0,dim(df)[1])
  v <- rep(0,dim(df)[1]) ;v.learn <- rep(0,dim(df)[1])
  z <- rep(0,dim(df)[1]) ;z.learn <- rep(0,dim(df)[1])
  n <- 0

  for(i in 1:(dim(df)[1]-1)){
    if(all(df[i,pa.names]==pa.val)){
      counts[as.numeric(df[i,node.idx])] <- counts[as.numeric(df[i,node.idx])] + 1
      new.alpha[as.numeric(df[i,node.idx])] <- new.alpha[as.numeric(df[i,node.idx])] + 1
      which.val <- as.numeric(df[(i+1),node.idx])
      p <- (alpha.vec)/sum(alpha.vec)
      n <- n+1

      s[i] <- -log(p)[which.val]
      e[i] <- -sum((p*log(p)))

      p.learn <- (new.alpha) /sum(new.alpha)
      #p.learn <- lgamma(sum(new.alpha))-sum(lgamma(new.alpha)) + sum(lgamma(new.alpha + counts)) - lgamma(sum(new.alpha + counts))
      s.learn[i] <- -log(p.learn)[which.val]

      e.learn[i] <- -sum((p.learn*log(p.learn)))
      pi <- e[i]/n; pi.learn <- e.learn[i]/n;

      v[i] <- (sum(p*log(p)^2)) -e[i]^2#+(pi^2/12)
      #v[i] <- n*pi*(1-pi)
      z[i] <- (cumsum(s)[i]-cumsum(e)[i])/sqrt(cumsum(v)[i])
      v.learn[i] <- (sum(p.learn*log(p.learn)^2)) -e.learn[i]^2#-(c^2/12)
      z.learn[i] <- (cumsum(s.learn)[i]-cumsum(e.learn)[i])/sqrt(cumsum(v.learn)[i])
    } else {
      new.alpha[as.numeric(df[i,node.idx])] <- new.alpha[as.numeric(df[i,node.idx])] + 1
      next
    }


  }
  #par(mfrow=c(1,2))
  #plot(s.learn)
  z.learn[which(z.learn==0)] <- NA
  return(z.learn)
}

#'@export
#' @importFrom ggplot2 ggtitle xlab ylab theme_minimal theme scale_colour_discrete
#' @importFrom reshape2 melt
#' @rdname sequential
pa.ch.graph <-  function(dframe, dag, node.name, pa.names, pa.val,alpha){
  z.learn <- sqntl.pa.ch.monitor(dframe,dag,node.name,pa.names,pa.val,alpha)
  t <- 1:length(z.learn[which(z.learn!=0)])
  pa.title <- toString(paste0(pa.names," = ",pa.val))
  data <- as.data.frame(cbind(t,z.learn[which(z.learn!=0)]))
  data <- melt(data,id='t')
  value <- data$value
  variable <- data$variable
  ggplot(data, aes(x=t, y=value, colour=variable)) + geom_point() + xlab('Relevant sample size') + ylab('Standardized Z Statistic') + theme_minimal() + ggtitle(paste0("p(",node.name," | ", pa.names, " =", pa.val, " ) ")) +   scale_colour_discrete(name="", labels=c( "Reference prior with learning")) + theme(legend.position="bottom")
}
