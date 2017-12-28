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
  \item{g.list}{A list of graph of class \code{igraph} with labelled vertices.}
}
\details{
  This function deals with bla bla bla
}
\value{
Returns a list of \code{dynsbm} objects. Each object of class \code{dynsbm} is a \code{list} with the following components:
\item{transition}{The Markov chain transition \code{matrix} of size (Q x Q).}
\item{membership}{An object of class \code{matrix} of size (N x T)
  containing the group membership estimated by MAP (>0, =0 for absent nodes).}
\item{beta}{An object of class \code{matrix} of size (T x Q x Q) containing
  the sparsity parameters of the model.}
\item{gamma}{Only if \code{edge.type=="discrete"}. An object of class \code{matrix} of size (T x Q x Q x K) containing
  the model parameters.}
\item{mu/sigma}{Only \code{edge.type=="continuous"}An object of class \code{matrix} of size (T x Q x Q)
  and a vector containing the model parameters.}
\item{loglikelihood}{Completed data log-likelihood.}
\item{iter}{Number of used algorithm iterations.}
\item{directed}{Specifies whether the model is build for directed networks.}
\item{self.loop}{Specifies whether the model allows self-loops.}
}
\references{
Catherine Matias and Vincent Miele, \emph{Statistical clustering of temporal networks through a dynamic stochastic block model}, Journal of the Royal Statistical Society: Series B (2016)
http://dx.doi.org/10.1111/rssb.12200
http://arxiv.org/abs/1506.07464

Vincent Miele and Catherine Matias, \emph{Revealing the hidden structure of dynamic ecological networks}, preprint
https://arxiv.org/abs/1701.01355
}
\author{
Authors: Catherine Matias, Vincent Miele
  
Maintainer: Vincent Miele <vincent.miele@univ-lyon1.fr>
}
\examples{
data(simdataT5Q4N40binary)

## estimation for Q=1..6 groups
list.dynsbm <- select.dynsbm(simdataT5Q4N40binary, 
				Qmin=1, Qmax=6, edge.type="binary", nstart=1)
				
\dontrun{
## better to use nstart>1 starting points
## but estimation can take 1-2 minutes
list.dynsbm <- select.dynsbm(simdataT5Q4N40binary, 
				Qmin=1, Qmax=6, edge.type="binary", nstart=25)
}

## selection of Q=4
dynsbm <- list.dynsbm[[4]]

## plotting intra/inter connectivity patterns
connectivity.plot(dynsbm, simdataT5Q4N40binary)

## plotting switches between groups
alluvial.plot(dynsbm)
}