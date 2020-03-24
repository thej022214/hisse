\name{makeGeoHiSSELikelihood}
\alias{makeGeoHiSSELikelihood}
\title{Export likelihood function for the GeoHiSSE model}
\description{Exports a likelihood function conditioned on the data and a
  named vector with the parameters for each of the models.}
\usage{
makeGeoHiSSELikelihood(phy, data, hidden.areas=0, f=c(1,1,1), assume.cladogenetic=TRUE, 
condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, dt.threads=1, ode.eps = 0, 
bad.likelihood = exp(-500))
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a data matrix containing species information (see Details).}
\item{hidden.areas}{Number of hidden areas included in the model. When hidden.areas = 0, the model is equal to the GeoSSE model (i.e., no hidden states included).}
\item{f}{vector of length 3 with the estimated proportion of extant species in areas 0, 1, and 2 (see fgeohisse function for more information) that are included in the phylogeny. By default all species are assumed to be sampled.}
\item{assume.cladogenetic}{a logical indicating whether the likelihood should be conditioned on geographical speciation happening at the nodes. If set to FALSE then all range transitions are treated as anagenetic. The default is \code{TRUE}. See Caetano et al. (2018 - Evolution) for more information.}
\item{condition.on.survival}{a logical indicating whether the likelihood should be conditioned on the survival of two lineages and the speciation event subtending them (Nee et al. 1994). The default is \code{TRUE}.}
\item{root.type}{indicates whether root summarization follow the procedure described by FitzJohn et al. 2009, \dQuote{madfitz} or Herrera-Alsina et al. 2018, \dQuote{herr_als}.}
\item{root.p}{a vector indicating fixed root state probabilities. The default is \code{NULL}.}
\item{dt.threads}{Number of core threads to be used by the data.table function. Default is 1.}
\item{ode.eps}{sets the tolerance for the integration at the end of a branch. Essentially if the sum of compD is less than this tolerance, then it assumes the results are unstable and discards them. The default is set to zero, but in testing a value of 1e-8 can sometimes produce stable solutions for both easy and very difficult optimization problems.}
\item{bad.likelihood}{Value returned when there is a problem in the computation of the likelihood value for the function.}
}
\details{
This function sets up and returns the likelihood for the GeoHiSSE model
together with a vector of parameters. The likelihood function is
conditioned on the observed data and will return a value of
loglikelihood given a vector of parameter values. The length of the
parameter vector as well as the order of
the parameter vector cannot be changed. Please pay special attention to the length of the parameter
vector and the names of the parameters provided by the \dQuote{pars}
element of the list.

Note that when the likelihood computation fails the function will return the value set as \dQuote{bad.likelihood}.
}
\value{
\code{makeGeoHiSSELikelihood} returns a list with elements:
\item{$loglik}{the likelihood function for the model. This has a single
  parameter 'p'.}
\item{$pars}{a named vector for the likelihood function pupulated with 0
  values.}
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
Caetano, D.S., B.C. O'Meara, and J.M. Beaulieu. 2018. Hidden state models improve state-dependent diversification approaches, including biogeographic models. Evolution, 72:2308-2324.

Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 

FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Syst. Biol. 58:595-611.

Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary characters effect on speciation and extinction. Syst. Biol. 56:701-710.

Nee S., May R.M., and Harvey P.H. 1994. The reconstructed evolutionary process. Philos. Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
}

\author{Jeremy M. Beaulieu and Daniel S. Caetano}
\keyword{models}
