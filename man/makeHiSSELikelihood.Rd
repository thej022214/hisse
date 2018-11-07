\name{makeHiSSELikelihood}
\alias{makeHiSSELikelihood}
\title{Export likelihood function for the HiSSE model}
\description{Exports a likelihood function dependent on the data and a
  named vector with the parameters for the model.}
\usage{
makeHiSSELikelihood(phy, data, f=c(1,1), condition.on.survival=TRUE, root.type="madfitz", 
root.p=NULL, ode.eps=0)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a data matrix containing species information (see Details).}
\item{f}{vector of length 2 with the estimated proportion of extant species in state 0 and 1 that are included in the phylogeny. A value of c(0.25, 0.5) means that 25 percent of species in state 0 and 50 percent of species in state 1 are included in the phylogeny. By default all species are assumed to be sampled.}
\item{condition.on.survival}{a logical indicating whether the likelihood should be conditioned on the survival of two lineages and the speciation event subtending them (Nee et al. 1994). The default is \code{TRUE}.}
\item{root.type}{indicates whether root summarization follow the procedure described by FitzJohn et al. 2009, \dQuote{madfitz} or Herrera-Alsina et al. 2018, \dQuote{herr_als}.}
\item{root.p}{a vector indicating fixed root state probabilities. The default is \code{NULL}.}
\item{ode.eps}{sets the tolerance for the integration at the end of a branch. Essentially if the sum of compD is less than this tolerance, then it assumes the results are unstable and discards them. The default is set to zero, but in testing a value of 1e-8 can sometimes produce stable solutions for both easy and very difficult optimization problems.}
}
\details{
This function sets up and returns the likelihood for the HiSSE model
together with a vector of parameters. The likelihood function is
conditioned on the observed data and will return a value of
loglikelihood given a vector of parameter values. The vector of
parameter values is populated with a suggested starting point. This is
the same starting point used by the 'hisse' function if none is provided
by the user. The length of the parameter vector as well as the order of
the parameter vector cannot be changed. The values of the parameter
vector need to be log-transformed.
}
\value{
\code{makeHiSSELikelihood} returns a list with elements:
\item{$loglik}{the likelihood function for the model. This has a single
  parameter 'p'.}
\item{$pars}{the named vector for the likelihood function. Values are log-transformed.}
}
\examples{
\donttest{
library(diversitree)
library(hisse)
## Generate data:
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(4) 
phy <- tree.bisse(pars, max.t=30, x0=0)
sim.dat <- data.frame(names(phy$tip.state), phy$tip.state)
## Get lik function:
model <- makeHiSSELikelihood(phy = phy, data = sim.dat)
names( model ) ## Components of the model object.
## Note that the parameters are log-transformed.
exp( model$pars )
model$loglik(p = model$pars) ## Compute the log likelihood.
}
}
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 

FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Syst. Biol. 58:595-611.

Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary characters effect on speciation and extinction. Syst. Biol. 56:701-710.

Nee S., May R.M., and Harvey P.H. 1994. The reconstructed evolutionary process. Philos. Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
}

\author{Jeremy M. Beaulieu}
\keyword{models}
