\name{hisse.new}
\alias{hisse.new}
\alias{fhisse}
\title{Hidden State Speciation and Extinction}
\description{Sets up and executes a HiSSE model (Hidden State Speciation and Extinction) 
on a phylogeny and character distribution.}
\usage{
hisse.new(phy, data, f=c(1,1), turnover=c(1,2), eps=c(1,2), 
hidden.states=FALSE, trans.rate=NULL, condition.on.survival=TRUE, 
root.type="madfitz", root.p=NULL, sann=FALSE, sann.its=10000, 
bounded.search=TRUE,  max.tol=.Machine$double.eps^.50, starting.vals=NULL,
turnover.upper=10000, eps.upper=3, trans.upper=100, restart.obj=NULL, 
ode.eps=0, dt.threads=1)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and
  with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a matrix (or dataframe) with three columns. The first column 
  containing the species names and the second and third containing the binary 
  character information. Character "0" is on column 2 and chracter "1" is on 
  column 3. A value of 0 means character absent and a value of 1 character 
  present. The input of data follows a Pagel model. See 'Details'.}
\item{f}{vector of length 2 with the estimated proportion of extant species in 
  state 0 and 1 that are included in the phylogeny. A value of c(0.25, 0.5) means 
  that 25 percent of species in state 0 and 50 percent of species in state 1 are 
  included in the phylogeny. By default all species are assumed to be sampled.}
\item{turnover}{a numeric vector indicating the number of free turnover parameters in the
  model.}
\item{eps}{a numeric vector indicating the number of free extinction fraction parameters 
  in the model.}
\item{hidden.states}{a logical indicating whether the model includes a
  hidden states. The default is \code{FALSE}.} 
\item{trans.rate}{provides the transition rate model. See function
  \code{TransMatMakerMuHiSSE}.}
\item{condition.on.survival}{a logical indicating whether the likelihood
  should be conditioned on the survival of two lineages and the
  speciation event subtending them (Nee et al. 1994). The default is \code{TRUE}.}
\item{root.type}{indicates whether root summarization follow the procedure described by 
  FitzJohn et al. 2009, \dQuote{madfitz} or Herrera-Alsina et al. 2018, \dQuote{herr_als}.}
\item{root.p}{a vector indicating fixed root state probabilities. The
  default is \code{NULL}.}
\item{sann}{a logical indicating whether a two-step optimization
  procedure is to be used. The first includes a simulate annealing
  approach, with the second involving a refinement using
  \code{subplex}. The default is \code{FALSE}.}
\item{sann.its}{a numeric indicating the number of times the simulated
  annealing algorithm should call the objective function.}
\item{bounded.search}{a logical indicating whether or not bounds should
  be enforced during optimization. The default is is \code{TRUE}.}
\item{max.tol}{supplies the relative optimization tolerance to
  \code{subplex}.}
\item{starting.vals}{a numeric vector of length 3 with starting values for the
  model for all areas and hidden states. Position [1] sets turnover, [2]
  sets extinction fraction, and [3] transition rates.}
\item{turnover.upper}{sets the upper bound for the turnover parameters.}
\item{eps.upper}{sets the upper bound for the eps parameters.}
\item{trans.upper}{sets the upper bound for the transition rate parameters.}
\item{restart.obj}{an object of class  that contains everything to restart an optimization.}
\item{ode.eps}{sets the tolerance for the integration at the end of a
  branch. Essentially if the sum of compD is less than this tolerance,
  then it assumes the results are unstable and discards them. The
  default is set to zero, but in testing a value of 1e-8 can sometimes
  produce stable solutions for both easy and very difficult optimization problems.}
\item{dt.threads}{sets the number of threads available to data.table. In practice this 
  need not change from the default of 1 thread, as we have not seen any speedup from 
  allowing more threads.}
}
\details{
This function sets up and executes a faster version of the HiSSE model. Note that the
four-state character-independent model can be called from this command in addition to the 
two-state BiSSE model and the full character-dependent HiSSE model. For details on the
various options see \code{hisse}. 

Also note that this code will replace the original hisse function in the next version.
}
\value{
\code{hisse.new} returns an object of class \code{fhisse.fit}. This is a list with 
elements:
\item{$loglik}{the maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample-size.}
\item{$solution}{a matrix containing the maximum likelihood estimates of the model 
parameters.}
\item{$index.par}{an index matrix of the parameters being estimated.}
\item{$f}{user-supplied sampling frequencies.}
\item{$hidden.states}{a logical indicating whether hidden states were included in the 
model.}
\item{$condition.on.surivival}{a logical indicating whether the
  likelihood was conditioned on the survival of two lineages and the speciation event 
subtending them.}
\item{$root.type}{indicates the user-specified root prior assumption.}
\item{$root.p}{indicates whether the user-specified fixed root probabilities.}
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
\donttest{
library(diversitree)
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(4) 
phy <- tree.bisse(pars, max.t=30, x0=0)
sim.dat <- data.frame(names(phy$tip.state), phy$tip.state)

## Fit BiSSE equivalent:
trans.rates.bisse <-  TransMatMakerfHiSSE(hidden.traits=0)
pp.bisse <- hisse.new(phy, sim.dat, hidden.states=FALSE, turnover=c(1,2), 
eps=c(1,2), trans.rate=trans.rates.bisse)

## Now fit HiSSE equivalent with a hidden state for state 1:
trans.rates.hisse <- TransMatMakerfHiSSE(hidden.traits=1)
pp.hisse <- hisse.new(phy, sim.dat, hidden.states=TRUE, turnover=c(1,2,1,2), 
eps=c(1,2,1,2), trans.rate=trans.rates.hisse)
}
}
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models 
of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 

FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating trait-dependent speciation 
and extinction rates from incompletely resolved phylogenies. Syst. Biol. 58:595-611.

Herrera-Alsina, L., P. van Els, and R.S. Etienne. 2018. Detecting the dependence of diversification on multiples traits from phylogenetic trees and trait data. Systematic Biology, 68:317-328.

Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary characters effect on 
speciation and extinction. Syst. Biol. 56:701-710.

Nee S., May R.M., and Harvey P.H. 1994. The reconstructed evolutionary process. Philos. 
Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
}

\author{Jeremy M. Beaulieu}
\keyword{models}
