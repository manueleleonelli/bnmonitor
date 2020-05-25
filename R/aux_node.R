standardize <- function(vec,j){#worst level is the actual observation, we're looking at the 'surprise' of seeing that observation.
  p <- unlist(vec) #pick the worst one
  sm <- -log(p)
  em <- - sum(p*log(p))
  vm <- sum(p*(log(p)^2)) - em^2
  zm <- sm[j]-em/sqrt(vm)
  return(zm)
}

pass.ev <-function(i, df, dag.grain){
  ev <- querygrain(setEvidence(dag.grain,evidence=list(df[dim(df)[1],-i]), nodes=colnames(df)[i]))
  ev2 <- ev[match(names(df),names(ev))]  #reordering to match
  evi <- ev2[[i]]
  return(evi)
}
