\name{getMetaweb}
\alias{getMetaweb}
\title{
Get metaweb 
}
\description{
Computation of the binary metaweb from a list of graph
}
\usage{
getMetaweb(graphList)
}
\arguments{
  \item{graphList}{A \code{list} of graph objects of class \code{igraph}.}
}
\details{
  This function computes the metaweb from a list of graph. It computes the union (in the sense of graph theory) of the set of graphs.
}
\value{
\code{getMetaweb} returns an object of class \code{igraph}.
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
## generating a list of igraph objects
graphList <- c()
for(i in 1:4){
  graphLocal <- erdos.renyi.game(60, type = 'gnp', p.or.m=0.1, directed=T)
  V(graphLocal)$name <- as.character(1:60)
  graphList <- c(graphList, list(graphLocal))
}
names(graphList) <- c("A","B","C","D")

## building the metaweb
graphMetaweb <- getMetaweb(graphList)
}