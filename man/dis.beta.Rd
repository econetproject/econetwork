\name{dis.beta}
\alias{dis.beta}
\title{
Computation of the dissimilarity matrix (pairwise beta-diversity) for a set of network
}
\description{
  Computation of the dissimilarity matrix for a set of network. Each value of the matrix is the pairwise beta-diversity, computed using Hill numbers (Renyi entropy).
}
\usage{
dis.beta(graph.list, groups, eta=2, framework=c('RLC','Tu'), type=c('P','L','Pi'))
}
\arguments{
  \item{graph.list}{A \code{list} of graph objects of class
    \code{igraph}.}
  \item{groups}{A named vector of integers (from \code{1:k}) indicating the group to
    which each node belongs to. The length of \code{groups} must
    correspond to the number of different nodes present in
    \code{graph.list}. The names \code{names(groups)} must
    correspond to the nodes names in \code{graph.list}.}
  \item{eta}{eta}
  \item{framework}{Bla bla bla. This should be one of 'RLC' or 'Tu'}
  \item{type}{Bla bla bla. This should be one of 'P','L' or'Pi'}
}
\details{
  This function 
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
}
\author{
  Author: Stéphane Dray, Vincent Miele, Marc Ohlmann, Wilfried Thuiller
  
  Maintainer: Wilfried Thuiller <wilfried.thuiller@univ-grenoble-alpes.fr>
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