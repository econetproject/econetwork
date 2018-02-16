\name{get.metaweb}
\alias{get.metaweb}
\title{
Get metaweb 
}
\description{
Computation of the binary metaweb from a list of graph
}
\usage{
get.metaweb(graph.list)
}
\arguments{
  \item{graph.list}{A \code{list} of graph objects of class \code{igraph}.}
}
\details{
  This function computes the metaweb from a list of graph. It computes the union (in the sense of graph theory) of the set of graphs.
}
\value{
\code{get.metaweb} returns an object of class \code{igraph}.
}
\references{

}
\author{
Authors: Stéphane Dray, Vincent Miele, Marc Ohlmann, Wilfried Thuiller
  
Maintainer: Wilfried Thuiller <wilfried.thuiller@univ-grenoble-alpes.fr>
}
\examples{
## generating a list of igraph objects
graph.list <- c()
for(i in 1:4){
  graph.local <- erdos.renyi.game(60, type = 'gnp', p.or.m=0.1, directed=T)
  V(graph.local)$name <- as.character(1:60)
  graph.list <- c(graph.list, list(graph.local))
}
names(graph.list) <- c("A","B","C","D")

## building the metaweb
graph.metaweb <- get.metaweb(graph.list)
}