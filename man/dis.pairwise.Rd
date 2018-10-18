\name{dis.pairwise}
\alias{dis.pairwise}
\title{
Computation of the dissimilarity matrix (pairwise beta-diversity) for a set of network
}
\description{
  Computation of the dissimilarity matrix for a set of network. Each value of the matrix is the pairwise beta-diversity, computed using Hill numbers. It measures either group, links, or probability of links dissimilarity.  
}
\usage{
dis.beta(graph.list, groups, eta=1, framework=c('RLC','Tu'), type=c('P','L','Pi'))
}
\arguments{
  \item{graph.list}{A \code{list} of graph objects of class
    \code{igraph}.}
  \item{groups}{A named vector of class \code{character}  indicating the group to
    which each node belongs to. The length of \code{groups} must correspond to the number of different nodes present in \code{graph.list}. The names \code{names(groups)} must
    correspond to the nodes names in \code{graph.list}.}
  \item{eta}{A positive number that controls the weight given to abundant groups/links. Default value is 1.}
  \item{type}{The type of diversity used to measure dissimilarity. It can be groups diversity ('P'), link diversity ('L') or probability of link diversity ('Pi').}
}
\details{
  
}
\value{
  Return a \code{matrix} whose elements are the pairwise dissimilarities. 
}
\references{

}
\author{
  Author: Stéphane Dray, Vincent Miele, Marc Ohlmann, Wilfried Thuiller
  Maintainer: Wilfried Thuiller <wilfried.thuiller@univ-grenoble-alpes.fr>
}
\examples{
#Generating a set of Erdos-Renyi graphs and give name to nodes.
n.graph=10
graph.list=c()
n=57 #number of nodes of each graph
C=0.1  #connectance of each graph
for(i in 1:n.graph){
  graph.loc=erdos.renyi.game(n,type = 'gnp',p.or.m = C,directed = T)
  V(graph.loc)$name=as.character(1:57)
  graph.list=c(graph.list,list(graph.loc))
}

g.list=graph.list
groups=c(rep("a",23),rep("b",34)) #vector that gives the group of each node
names(groups)=as.character(1:57)

dis.pairwise(graph.list,groups,type = 'L') #dissimilarity matrix based on links beta-diversity
}