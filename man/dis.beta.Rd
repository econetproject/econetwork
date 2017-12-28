\name{dis.beta}
\alias{dis.beta}
\title{
Computation of the dissimilairty matrix (pairwise beta-diversity) for a set of network
}
\description{
Computation of the dissimilarity matrix for a set of network. Each value of the matrix is the pairwise beta-diversity, computed using Hill numbers (Renyi entropy).
}
\usage{
dis.beta(g.list, groups, eta=2, framework=c('RLC','Tu'), type=c('P','L','Pi'))
}
\arguments{
  \item{g.list}{A list of graph of class \code{igraph} with labelled vertices.}
  \item{groups}{An object of class \code{matrix} of size (N x 2), where  N is the number of species in the metaweb. The first column of this matrix contains the names
  of the species present in the metaweb whereas the second column contains the group of the species present in the metaweb. The default value is a matrix that contains two identical columns :
      the names of the species in the metaweb}
  \item{eta}{The order of the Hill number (also knwon as Renyi entropy) used to compute the beta-diversity}
  \item{Qmax}{Maximum number of groups.}
  \item{edge.type}{Type of adjacency matrices. This should be (an unambiguous abbreviation of) one of \code{binary}, \code{discrete} or \code{continuous}. 
  See the "Details" section.}
  \item{K}{Only if \code{edge.type=="discrete"}. Number of non-zero discrete values (i.e. in {1,..,K}).} 
  \item{directed}{Whether to consider the network is directed (and
    therefore whether \code{Y} is asymmetric).}
  \item{self.loop}{Whether to allow and consider self-loops (edges from one node
    to the same node) in the estimation steps.}   
  \item{nb.cores}{Number of cores to use, i.e. how many child
    processes and how many threads will be run simultaneously during the
    initialization and the estimation steps respectively.}
  \item{iter.max}{Maximal number of algorithm iterations.}
  \item{nstart}{Number of starting points for the iterative estimation algorithm.  See the "Details" section.}
  \item{perturbation.rate}{Rate of perturbation (in [0,1], see \code{nstart}) for the iterative estimation algorithm. This rate is
    the fraction of nodes for which its group is randomly shuffled.}
  \item{plot}{Display a plot with the loglikelihood and the ICL criteria if \code{edge.type=="binary"/"continuous"}. See the "Details" section.} 
}
\details{
  This function deals with binary or weighted dynamic/temporal/evolving networks (with discrete or continuous edges).
  The adjacency matrices must be coded with 0/1 in the \code{binary} case, 
  with 0/y where y belongs to the set {1,..,K} in the \code{discrete} case 
  or with 0/y where y is numeric, must be positive and is supposed to fit a gaussian mixture in
  the \code{continuous} case.

  Presence/absence information allows to model node's arrival or
  departure, birth or death, or simply enables to specify missing data
  (as absent nodes). If this information is missing (NULL), the presence/absence is deduced automatically from \code{Y} 
    by searching for nodes that do not participate in any edges (lines/columns of O in \code{Y}) and declaring them as absent. 
  This function does not support the existence of nodes that are never present (error message in this case).  
  
  The estimation algorithm is iterative and rely on a starting point. Therefore, it is possible to start the algorithm many times with 'nstart' starting points. 
  The first starting point is obtained with an ad-hoc use of the \code{kmeans} function. 
  The follwing starting point are obtained by perturbating the first one (see \code{perturbation.rate}).
  The greater \code{nstart}, the more accurate the results.
  
  To select the best number of groups, the "elbow" method consists in finding the point where the slope of the loglikelihood is decreasing (i.e. the loglikelihood is reaching a plateau).
  If \code{edge.type=="binary"}, the ICL criteria (plotted in red) has to be used : the best number of nodes is supposed to maximize the ICL criteria. 
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