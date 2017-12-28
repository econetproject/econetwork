\name{get.metaweb}
\alias{get.metaweb}
\title{
Get metaweb bla bla bla
}
\description{
Computation of bla bla bla
}
\usage{
get.metaweb(graph.list)
}
\arguments{
  \item{graph.list}{A \code{list} of graph objects of class \code{igraph}.}
}
\details{
  This function  bla bla bla
}
\value{
\code{get.metaweb} returns an object of class \code{igraph}.
}
\references{
Catherine Matias and Vincent Miele, \emph{Statistical clustering of temporal networks through a dynamic stochastic block model}, Journal of the Royal Statistical Society: Series B (2016)
http://dx.doi.org/10.1111/rssb.12200
http://arxiv.org/abs/1506.07464

Poisot et al. 2011
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