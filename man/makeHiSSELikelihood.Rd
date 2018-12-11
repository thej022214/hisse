\name{makeHiSSELikelihood}
\alias{makeHiSSELikelihood}
\title{Export likelihood function for the HiSSE model}
\description{Exports a likelihood function conditioned on the data and a
  named vector with the parameters for each of the models.}
\usage{
makeHiSSELikelihood(phy, data, hidden.states=TRUE, null4=FALSE,
f=c(1,1), condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, ode.eps=0)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a data matrix containing species information (see Details).}
\item{hidden.states}{whether the model has hidden states. This is
  ignored if the option \dQuote{null4} is selected.}
\item{null4}{used to select the null model for the full HiSSE
  model. This is a special model with 4 hidden states.}
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
loglikelihood given a vector of parameter values. The length of the
parameter vector as well as the order of
the parameter vector cannot be changed. The parameter values are
provided in natural form but are log-transformed for the likelihood
evaluation. Please pay special attention to the length of the parameter
vector and the names of the parameters provided by the \dQuote{pars}
element of the list.

When the option \dQuote{null4} is set to TRUE, then the likelihood
returned is for the HiSSE null model with 4 hidden states. The returned
list will include an additional element named \dQuote{trans.mat.guide}
which can be used as a reference for the meaning of the 32 transition
parameters. Note that the original model in \dQuote{hisse.null4} usually
sets all these transitions to a single value. This helps model
estimation since information is limited to estimate distinct transition
rates between parameters.
}
\value{
\code{makeHiSSELikelihood} returns a list with elements:
\item{$loglik}{the likelihood function for the model. This has a single
  parameter 'p'.}
\item{$pars}{a named vector for the likelihood function pupulated with 0
  values.}
\item{$trans.mat.guide}{a reference matrix to understand the parameters
  returned by the \dQuote{$pars} vector. Only present if the
  \dQuote{null4} option is set to TRUE.}
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
lik.hisse <- makeHiSSELikelihood(phy = phy, data = sim.dat, hidden.states = FALSE)
likf <- lik.hisse$log.lik
pars <- lik.hisse$pars
## Set the parameter values. Note that we have turnover and eps, not speciation and extinction!
pars <- setNames(c(0.1+0.03,0.2+0.03,0.03/0.1,0.03/0.2,0.01,0.01), names(pars))
## Compute the log-likelihood for the model.
likf(pars)
}
}
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 

FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Syst. Biol. 58:595-611.

Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary characters effect on speciation and extinction. Syst. Biol. 56:701-710.

Nee S., May R.M., and Harvey P.H. 1994. The reconstructed evolutionary process. Philos. Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
}

\author{Jeremy M. Beaulieu and Daniel S. Caetano}
\keyword{models}
