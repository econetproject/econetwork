\name{disPairwise}
\alias{disPairwise}
\title{
Computation of the dissimilarity matrix (pairwise beta-diversity) for a set of network
}
\description{
  Computation of the dissimilarity matrix for a set of network. Each value of the matrix is the pairwise beta-diversity, computed using Hill numbers. It measures either group, links, or probability of links dissimilarity.  
}
\usage{
disPairwise(graphList, groups=NULL, eta=1, type=c('P','L','Pi'))
}
\arguments{
  \item{graphList}{A \code{list} of graph objects of class
    \code{igraph}.}
  \item{groups}{A named vector of class \code{character}  indicating the group to
    which each node belongs to. The length of \code{groups} must correspond to the number of different nodes present in \code{graphList}. The names \code{names(groups)} must
    correspond to the nodes names in \code{graphList}. If NULL, the groups are the initial nodes.}
  \item{eta}{A positive number that controls the weight given to abundant groups/links. Default value is 1.}
  \item{type}{The type of diversity used to measure dissimilarity. It can be groups diversity ('P'), link diversity ('L') or probability of link diversity ('Pi').}
}
\details{
  
}
\value{
  Return a \code{matrix} whose elements are the pairwise dissimilarities. 
}
\references{
Marc Ohlmann, Vincent Miele, Stéphane Dray, Loïc Chalmandrier, Louise
O’Connor & Wilfried Thuiller, Diversity indices for ecological networks
: a unifying framework using Hill numbers. Ecology Letters (2019)
}
\author{
Authors: Stéphane Dray, Vincent Miele, Marc Ohlmann, Wilfried Thuiller
Maintainer: Wilfried Thuiller <wilfried.thuiller@univ-grenoble-alpes.fr>
}
\examples{
#Generating a set of Erdos-Renyi graphs and give name to nodes.
nbGraph=10
graphList=c()
n=57 #number of nodes of each graph
C=0.1  #connectance of each graph
for(i in 1:nbGraph){
  graphLoc = erdos.renyi.game(n, type = 'gnp', p.or.m = C, directed = T)
  V(graphLoc)$name = as.character(1:57)
  graphList = c(graphList,list(graphLoc))
}

gList=graphList
groups=c(rep("a",23),rep("b",34)) #vector that gives the group of each node
names(groups)=as.character(1:57)

disPairwise(graphList,groups,type = 'L') #dissimilarity matrix based on links beta-diversity
}