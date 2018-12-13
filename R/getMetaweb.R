get.metaweb <- function(g.list){#get the metaweb of a list of graph (as the union on a set of graph)
  g.metaweb <- g.list[[1]]
  for(i in 2:length(g.list)){
    g.metaweb <- g.metaweb %u% g.list[[i]]
  }
  return(g.metaweb)
}
