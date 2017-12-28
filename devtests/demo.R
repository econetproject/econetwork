library(igraph)
library(betalink)
library(divnet)

#We generate a set of Erdos-Renyi graphs and gives ids.
n.graph=10
graph.list=c()

set.seed(1)
n=57 #number of nodes of each graph, of course, we could have chosen a value for each graph
C=0.1  #connectance of each graph

groups <- rep(-1,n)
groups[1:19]=1
groups[20:38]=2
groups[39:57]=3
names(groups) <- 1:57 ## WARNING cette ligne est nécessaire

for(i in 1:n.graph){
  graph.loc=erdos.renyi.game(n,type = 'gnp',p.or.m = C,directed = T)
  V(graph.loc)$name=as.character(1:n)
  graph.list=c(graph.list,list(graph.loc))
}
names(graph.list)=LETTERS[1:10]

##source("../R/get.metaweb.R")
g.metaweb <- get.metaweb(graph.list)

##source("../R/utils.R")
meta.array <- metaweb.params(graph.list, groups)

##source("../R/div.partition.R")
div.partition(graph.list, groups, eta=1, framework = 'RLC',type = 'P')
div.partition(graph.list, groups, eta=1, framework = 'RLC',type = 'L')
div.partition(graph.list, groups, eta=1, framework = 'RLC',type = 'Pi')  #same value due to the fact that we have an equal proportion of ids...

##source("../R/chalmandrier.R")
div.partition(graph.list, groups, eta=1, framework = 'Tu',type = 'P')
div.partition(graph.list, groups, eta=1, framework = 'Tu',type = 'L')
div.partition(graph.list, groups, eta=1, framework = 'Tu',type = 'Pi')  #same value due to the fact that we have an equal proportion of ids...

#dissimilarity matrix

##source("../R/dis.beta.R")
dis.beta(graph.list, groups, eta=1, framework = 'RLC',type = 'P')
dis.beta(graph.list, groups, eta=1, framework = 'RLC',type = 'L')
dis.beta(graph.list, groups, eta=1, framework = 'RLC',type = 'Pi')

dis.beta(graph.list, groups, eta=1, framework = 'RLC',type = 'P')
dis.beta(graph.list, groups, eta=1, framework = 'RLC',type = 'L')
dis.beta(graph.list, groups, eta=1, framework = 'RLC',type = 'Pi')




