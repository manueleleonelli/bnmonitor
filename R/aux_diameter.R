tvd <- function(x,y) sum(abs(x-y))*0.5


diam <- function(bn, node){
  if(class(bn)[1] != "bn.fit"){stop("The input is not a bn.fit object")}
  if(!node%in% nodes(bn)){stop("Not a valid vertex name in input")}
  object <- bn[[node]]
  probi <- object$prob
  levels_out <- dim(probi)[1]
  tot <- length(probi)
  ind <- c(which(1:tot%% levels_out == 1),tot+1)
  max <- 0
  for(i in 1:(length(ind)-2)){
    for(j in (1+i):(length(ind)-1)){
      temp <- tvd(probi[ind[i]:(ind[i+1]-1)],probi[ind[j]:(ind[j+1]-1)])
      if(max < temp) max <- temp
    }
  }
  return(max)
}


mi <- function(bnfit, node1, node2){
  if(class(bnfit)[1] != "bn.fit"){stop("An input bn.fit object required")}
  if(!node1%in% nodes(bnfit)){stop("Not a valid vertex name in input")}
  if(!node2%in% nodes(bnfit)){stop("Not a valid vertex name in input")}
  try <- as.grain(bnfit)
  proby1 <- querygrain(try,node1)
  proby2 <- querygrain(try,node2)
  proby <- querygrain(try,c(node1,node2),type="joint")
  if(names(dimnames(proby))[1]==node1){probs <- outer(proby1[[1]],proby2[[1]])} else{
    probs <- outer(proby2[[1]],proby1[[1]])
  }
  return(sum(proby*log(proby) - proby*log(probs)))
}
