get.metaweb.R <- function(g.list){
  g.metaweb <- g.list[[1]]
  for(i in 2:length(gr.list)){
    g.metaweb <- g.metaweb %u% g.list[[i]]
  }
  return(g.metaweb)
}
