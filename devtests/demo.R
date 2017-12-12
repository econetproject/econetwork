library(igraph)
library(betalink)

#We generate a set of Erdos-Renyi graphs and gives ids.
n.graph=10
graph.list=c()

set.seed(1)
n=57 #number of nodes of each graph, of course, we could have chosen a value for each graph
C=0.1  #connectance of each graph

for(i in 1:n.graph){
  graph.loc=erdos.renyi.game(n,type = 'gnp',p.or.m = C,directed = T)
  V(graph.loc)$name=as.character(1:n)
  V(graph.loc)$id[1:19]=1
  V(graph.loc)$id[20:38]=2
  V(graph.loc)$id[39:57]=3
  graph.list=c(graph.list,list(graph.loc))
}
names(graph.list)=LETTERS[1:10]

source("../R/get.metaweb.R")
g.metaweb <- get.metaweb(graph.list)

#get the metaweb array
meta.array <- metaweb(graph.list)

#diversity index
ABG_decomp_eta(meta.array,eta = 1,framework = 'RLC',type = 'P')
ABG_decomp_eta(meta.array,eta = 1,framework = 'RLC',type = 'L')
ABG_decomp_eta(meta.array,eta = 1,framework = 'RLC',type = 'Pi')  #same value due to the fact that we have an equal proportion of ids...

ABG_decomp_eta(meta.array,eta = 1,framework = 'Tu',type = 'P')
ABG_decomp_eta(meta.array,eta = 1,framework = 'Tu',type = 'L')
ABG_decomp_eta(meta.array,eta = 1,framework = 'Tu',type = 'Pi')  #same value due to the fact that we have an equal proportion of ids...

#dissimilarity matrix

Beta_dis_eta(metaweb.array = meta.array,eta=1,framework = 'RLC',type = 'P')
Beta_dis_eta(metaweb.array = meta.array,eta=1,framework = 'RLC',type = 'L')
Beta_dis_eta(metaweb.array = meta.array,eta=1,framework = 'RLC',type = 'Pi')

Beta_dis_eta(metaweb.array = meta.array,eta=1,framework = 'RLC',type = 'P')
Beta_dis_eta(metaweb.array = meta.array,eta=1,framework = 'RLC',type = 'L')
Beta_dis_eta(metaweb.array = meta.array,eta=1,framework = 'RLC',type = 'Pi')




