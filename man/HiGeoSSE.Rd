\name{HiGeoSSE}
\alias{HiGeoSSE}
\title{Hidden Geographic State Speciation and Extinction}
\description{Sets up and executes a HiGeoSSE model (Hidden Geographic State Speciation and Extinction) on a phylogeny and character distribution.}
\usage{
HiGeoSSE(phy, data, f=c(1,1,1), speciation=c(1,2,3), extirpation=c(1,2), hidden.areas=FALSE, trans.rate=NULL, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, sann=FALSE, sann.its=10000, bounded.search=TRUE, max.tol=.Machine$double.eps^.25, mag.san.start=0.5, starting.vals=NULL, speciation.upper=1000, extirpation.upper=1000, trans.upper=100, ode.eps=0)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a data matrix containing species information (see Details).}
\item{f}{vector of length 3 with the estimated proportion of extant
  species in areas 0 (or AB), 1 (or A), and 2 (or B) that are included
  in the phylogeny. A value of c(0.25, 0.25, 0.5) means that 25 percent
  of species in areas AB and A and 50 percent of species in area B are included in the phylogeny. By default all species are assumed to be sampled.}
\item{hidden.areas}{a logical indicating whether the model includes a hidden area. The default is \code{FALSE}.} 
\item{trans.rate}{provides the transition rate model. See function \code{TransMatMakerHiGeoSSE}.}
\item{condition.on.survival}{a logical indicating whether the likelihood should be conditioned on the survival of two lineages and the speciation event subtending them (Nee et al. 1994). The default is \code{TRUE}.}
\item{root.type}{indicates whether root prior assumption should based the procedure described by FitzJohn et al. 2009, \dQuote{madfitz}, assumed equal, \dQuote{equal}, or set to user, \dQuote{user}.}
\item{root.p}{a vector indicating fixed root state probabilities. The default is \code{NULL}.}
\item{sann}{a logical indicating whether a two-step optimization procedure is to be used. The first includes a simulate annealing approach, with the second involving a refinement using \code{subplex}. The default is \code{FALSE}.}
\item{sann.its}{a numeric indicating the number of times the simulated annealing algorithm should call the objective function.}
\item{bounded.search}{a logical indicating whether or not bounds should be enforced during optimization. The default is is \code{TRUE}.}
\item{max.tol}{supplies the relative optimization tolerance to
  \code{subplex}.}
\item{mag.san.start}{UNCLEAR WHAT THIS DOES! ASK JEREMY!}
\item{starting.vals}{a vector of starting values to be used instead of
  the default settings. NEED MORE INFO ON THIS!}
\item{speciation.upper}{sets the upper bound for the speciation parameters.}
\item{extirpation.upper}{sets the upper bound for the extirpation parameters.}
\item{trans.upper}{sets the upper bound for the transition rate parameters.}
\item{ode.eps}{sets the tolerance for the integration at the end of a branch. Essentially if the sum of compD is less than this tolerance, then it assumes the results are unstable and discards them. The default is set to zero, but in testing a value of 1e-8 can sometimes produce stable solutions for both easy and very difficult optimization problems.}
}
\details{
This function sets up and executes the HiGeoSSE model. The model closely
follows \code{diversitree}, although here we employ modified
optimization procedures. As for data file format, \code{HiGeoSSE}
expects a two column matrix or data frame, with the first column
containing the species names and the second containing the are
information. The area information need to be in the format of three
numbers: 0 for widespread area (i.e., AB), 1 for area A, and 2 for area
B. Note that the order of the data file and the names in the
\dQuote{phylo} object need not be in the same order; \code{hisse} deals
with this internally. Also, the character information MUST be one of 0,
1 or 2, otherwise, the function will misbehave.

To setup a model, users input vectors containing values to indicate how
many free parameters are to be estimated for each of the variables in
the model. This is very similar to the procedure used to set a
\code{hisse} model. Once the model is specified, the parameters can be estimated using the subplex routine (default), or use a two-step process (i.e., sann=TRUE) that first employs a stochastic simulated annealing procedure, which is later refined using the subplex routine.

The \dQuote{trans.rate} input is the transition model and has an entirely different setup than speciation and extirpation rates. See \code{TransMatMakerHiGeoSSE} function for more details. 

For user-specified \dQuote{root.p}, you should specify the probability
for each area. If you are doing a hidden model, there will be six areas:
0A, 1A, 2A, 0B, 1B, 2B. So if you wanted to say the root had to be in
area 0 (widespread distribution), you would specify \dQuote{root.p =
  c(0.5, 0, 0, 0.5, 0, 0)}. In other words, the root has a 50\% chance
to be in one of the areas 0A or 0B.
}
\value{
\code{HiGeoSSE} returns an object of class \code{higeosse.fit}. This is a list with elements:
\item{$loglik}{the maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample-size.}
\item{$solution}{a matrix containing the maximum likelihood estimates of the model parameters.}
\item{$index.par}{an index matrix of the parameters being estimated.}
\item{$f}{user-supplied sampling frequencies.}
\item{$hidden.areas}{a logical indicating whether hidden areas were included in the model.}
\item{$condition.on.surivival}{a logical indicating whether the likelihood was conditioned on the survival of two lineages and the speciation event subtending them.}
\item{$root.type}{indicates the user-specified root prior assumption.}
\item{$root.p}{indicates whether the user-specified fixed root probabilities.}
\item{$timeslice}{indicates whether the user-specified timeslice that split the tree.}
\item{$phy}{user-supplied tree}
\item{$data}{user-supplied dataset}
\item{$trans.matrix}{the user-supplied transition matrix}
\item{$max.tol}{relative optimization tolerance.}
\item{$starting.vals}{The starting values for the optimization.}
\item{$upper.bounds}{the vector of upper limits to the optimization search.}
\item{$lower.bounds}{the vector of lower limits to the optimization
  search.}
\item{$ode.eps}{The ode.eps value used for the estimation.}
}
\examples{
## NEED TO POPULATE THIS!
}
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 

FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Syst. Biol. 58:595-611.

Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary characters effect on speciation and extinction. Syst. Biol. 56:701-710.

Nee S., May R.M., and Harvey P.H. 1994. The reconstructed evolutionary process. Philos. Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
}

\author{Jeremy M. Beaulieu}
\keyword{models}
