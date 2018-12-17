library(econetwork)
library(igraph)

#We generate a set of Erdos-Renyi graphs and gives ids.
nbGraph = 10
gList = c()

set.seed(1)
n = 57 #number of nodes of each graph, of course, we could have chosen a value for each graph
C = 0.1  #connectance of each graph

groups = rep(-1,n)
groups[1:19] = 1
groups[20:38] = 2
groups[39:57] = 3
names(groups) = 1:57 ## WARNING cette ligne est nécessaire

for(i in 1:nbGraph){
  graphLoc = erdos.renyi.game(n,type = 'gnp',p.or.m = C,directed = T)
  V(graphLoc)$name = as.character(1:n)
  gList = c(gList,list(graphLoc))
}
names(gList) = LETTERS[1:10]

##source("../R/get.metaweb.R")
gMetaweb <- getMetaweb(gList)

TEST <- FALSE
if(TEST){
    library(Matrix.utils)
    source("../R/utils.R")
    metaArray <- metawebParams(gList, groups)
}

divPartition(gList, groups, eta=1, framework='RLC', type='P')
divPartition(gList, groups, eta=1, framework='RLC', type='L')
divPartition(gList, groups, eta=1, framework='RLC', type='Pi')  #same value due to the fact that we have an equal proportion of ids...

divPartition(gList, groups, eta=1, framework='Chao', type='P')
divPartition(gList, groups, eta=1, framework='Chao', type='L')
divPartition(gList, groups, eta=1, framework='Chao', type='Pi')  #same value due to the fact that we have an equal proportion of ids...

disPairwise(gList, groups, eta=1, type='P')
disPairwise(gList, groups, eta=1, type='L')
disPairwise(gList, groups, eta=1, type='Pi')




