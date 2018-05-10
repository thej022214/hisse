\name{MuHiSSE}
\alias{MuHiSSE}
\title{Multicharacter Hidden State Speciation and Extinction}
\description{Sets up and executes a MuHiSSE model (Multicharacter Hidden State Speciation
 and Extinction) on a phylogeny and character distribution.}
\usage{
MuHiSSE(phy, data, f=c(1,1,1,1), turnover=c(1,2,3,4), eps=c(1,2,3,4), 
hidden.states=FALSE, trans.rate=NULL, condition.on.survival=TRUE, 
root.type="madfitz", root.p=NULL, sann=FALSE, sann.its=10000, 
bounded.search=TRUE,  max.tol=.Machine$double.eps^.50, starting.vals=NULL,
turnover.upper=10000, eps.upper=3, trans.upper=100, ode.eps=0)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and
    with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a matrix (or dataframe) with two columns containing species
  information. First column has the species names and second column has
  state codes. Values for the areas need to be 0, 1, or 2, where 0 is the 
  widespread area '01', 1 is endemic area '0' and 2 is endemic area '1'. See 'Details'.}
\item{f}{vector of length 3 with the estimated proportion of extant
  species in areas 0 (or '01'), 1 (or '0'), and 2 (or '1') that are included
  in the phylogeny. A value of c(0.25, 0.25, 0.5) means that 25 percent
  of species in areas '01' and '0' and 50 percent of species in area '1' are
  included in the phylogeny. By default all species are assumed to be
  sampled.}
\item{turnover}{a numeric vector of length equal to 4+(number of
  \code{hidden.states} * 4). A MuSSE model has 4 speciation parameters:
  lambda00, lambda01, lambda10, and lambda11. A MuHiSSE model with 
  one hidden states has 8 speciation parameters: lambda0A, s1A, s01A, s0B, s1B, s01B. And so on. The length of the numeric
vector needs to match the number of speciation parameters in the
model. See 'Details'.}
\item{eps}{a numeric vector of length equal to4+(number of
  \code{hidden.states} * 4). A MuSSE model has 4 extinction parameters:
  mu00, mu01, mu10, and mu11. A MuHiSSE model with one hidden area has 8 extinction
parameters: mu00A, mu01A, mu10A, mu11A, mu00B, mu01B, mu10B, and mu11B. And so on. The length of the numeric
vector needs to match the number of extinction parameters in the
model. See 'Details'.}
\item{hidden.states}{a logical indicating whether the model includes a
  hidden states. The default is \code{FALSE}.} 
\item{trans.rate}{provides the transition rate model. See function
  \code{TransMatMakerMuHiSSE}.}
\item{condition.on.survival}{a logical indicating whether the likelihood
  should be conditioned on the survival of two lineages and the
  speciation event subtending them (Nee et al. 1994). The default is \code{TRUE}.}
\item{root.type}{indicates whether root prior assumption should based
  the procedure described by FitzJohn et al. 2009, \dQuote{madfitz},
  assumed equal, \dQuote{equal}, or set to user, \dQuote{user}.}
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
\item{starting.vals}{a vector of starting values for the diversification parameters
  to be used instead of the default settings.}
\item{turnover.upper}{sets the upper bound for the turnover parameters.}
\item{eps.upper}{sets the upper bound for the eps parameters.}
\item{trans.upper}{sets the upper bound for the transition rate parameters.}
\item{ode.eps}{sets the tolerance for the integration at the end of a
  branch. Essentially if the sum of compD is less than this tolerance,
  then it assumes the results are unstable and discards them. The
  default is set to zero, but in testing a value of 1e-8 can sometimes
  produce stable solutions for both easy and very difficult optimization problems.}
}
\details{
MORE HERE.

The \dQuote{trans.rate} input is the transition model and has an
entirely different setup than speciation and extinction rates.
See \code{TransMatMakerMuHiSSE} function for more details. 

For user-specified \dQuote{root.p}, you should specify the probability
for each area. If you are doing a hidden model, there will be six areas:
0A, 1A, 2A, 0B, 1B, 2B. So if you wanted to say the root had to be in
area 0 (widespread distribution), you would specify \dQuote{root.p =
  c(0.5, 0, 0, 0.5, 0, 0)}. In other words, the root has a 50\% chance
to be in one of the areas 0A or 0B.
}
\value{
\code{MuHiSSE} returns an object of class \code{muhisse.fit}. This is a list with 
elements:
\item{$loglik}{the maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample-size.}
\item{$solution}{a matrix containing the maximum likelihood estimates of the model 
parameters.}
\item{$index.par}{an index matrix of the parameters being estimated.}
\item{$f}{user-supplied sampling frequencies.}
\item{$hidden.areas}{a logical indicating whether hidden areas were included in the 
model.}
\item{$condition.on.surivival}{a logical indicating whether the
  likelihood was conditioned on the survival of two lineages and the speciation event 
subtending them.}
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
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models 
of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 

FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating trait-dependent speciation 
and extinction rates from incompletely resolved phylogenies. Syst. Biol. 58:595-611.

Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary characters effect on 
speciation and extinction. Syst. Biol. 56:701-710.

Nee S., May R.M., and Harvey P.H. 1994. The reconstructed evolutionary process. Philos. 
Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
}

\author{Jeremy M. Beaulieu}
\keyword{models}
