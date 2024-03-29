\name{MiSSE}
\alias{MiSSE}
\title{Character-free State Speciation and Extinction}
\description{Sets up and executes a MiSSE model (Missing State Speciation
 and Extinction) on a phylogeny.}
\usage{
MiSSE(phy, f=1, turnover=c(1,2), eps=c(1,2), fixed.eps=NULL, condition.on.survival=TRUE,
root.type="madfitz", root.p=NULL, includes.fossils=FALSE, k.samples=NULL, 
strat.intervals=NULL, sann=TRUE, sann.its=5000, sann.temp=5230, sann.seed=-100377, 
bounded.search=TRUE, max.tol=.Machine$double.eps^.50, starting.vals=NULL, 
turnover.upper=10000, eps.upper=3, trans.upper=100, restart.obj=NULL, ode.eps=0, 
dt.threads=1, expand.mode=FALSE)
}
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format. If \code{includes.fossils=TRUE} then the input phy object must include extinct tips.}
\item{f}{the estimated proportion of extant species included in the phylogeny. A value of
  0.50 means that 50 percent of species are contained in the. By default all species are assumed
  to be sampled.}
\item{turnover}{a numeric vector of length equal to the number of suspected rates
  in turnover. See 'Details'.}
\item{eps}{a numeric vector of length equal to the number of suspected rates
  in extinction fraction. See 'Details'.}
\item{fixed.eps}{a value to be used to fix extinction fraction during search. Default is \code{NULL} meaning that it is freely estimated.}
\item{condition.on.survival}{a logical indicating whether the likelihood
  should be conditioned on the survival of two lineages and the
  speciation event subtending them (Nee et al. 1994). The default is \code{TRUE}.}
\item{root.type}{indicates whether root summarization follow the procedure described by FitzJohn et al. 2009, \dQuote{madfitz} or Herrera-Alsina et al. 2018, \dQuote{herr_als}.}
\item{root.p}{a vector indicating fixed root state probabilities. The
  default is \code{NULL}.}
\item{includes.fossils}{a logical indicating whether the tree contains fossil taxa. The default is \code{FALSE}.}
\item{k.samples}{a table of extinct individuals with sampled descendants. See vignette for how the table must be formatted.}
\item{strat.intervals}{a table of extinct individuals with sampled descendants. See vignette for how the table must be formatted.}
\item{sann}{a logical indicating whether a two-step optimization
  procedure is to be used. The first includes a simulate annealing
  approach, with the second involving a refinement using
  \code{subplex}. The default is \code{TRUE}.}
\item{sann.its}{a numeric indicating the number of times the simulated
  annealing algorithm should call the objective function.}
\item{sann.temp}{the starting temperature for the simulated annealing. Higher temperatures results in the chain sampling a much wider space initially. The default of 5320 is based on the default of the GenSA package. For larger trees setting this value higher in conjunction with more sann.its can drastically improve performance.}
\item{sann.seed}{the seed number for the simulated annealing algorithm. This value must be negative and an odd number.}
\item{bounded.search}{a logical indicating whether or not bounds should
  be enforced during optimization. The default is \code{TRUE}.}
\item{max.tol}{supplies the relative optimization tolerance to
  \code{subplex}.}
\item{starting.vals}{a numeric vector of length 3 with starting values for the
  model. Position [1] sets turnover, [2] sets extinction fraction, and
  [3] transition rates between distinct diversification rates.}
\item{turnover.upper}{sets the upper bound for the turnover parameters.}
\item{eps.upper}{sets the upper bound for the eps parameters.}
\item{trans.upper}{sets the upper bound for the transition rate parameters.}
\item{restart.obj}{an object of class  that contains everything to restart an optimization.}
\item{ode.eps}{sets the tolerance for the integration at the end of a
  branch. Essentially if the sum of compD is less than this tolerance,
  then it assumes the results are unstable and discards them. The
  default is set to zero, but in testing a value of 1e-8 can sometimes
  produce stable solutions for both easy and very difficult optimization problems.}
\item{dt.threads}{sets the number of threads available to data.table. In practice this need not change from the default of 1 thread, as we have not seen any speedup from allowing more threads.}
\item{expand.mode}{allows passing in the number of free parameters for turnover and eps and creates vectors to permit this.}
}
\details{
One thing pointed out in the original HiSSE paper (Beaulieu & O'Meara, 2016) is
that the trait-independent hisse model is basically a model for traits and a separate
model for shifts in diversification parameters, much like BAMM (though without priors,
discontinuous inheritance of extinction probability, or other mathematical foibles).
The hidden states can drive different diversification processes, and the traits just
evolve in a regular trait model. At that point, there is no harm in just dropping the
trait (or analyzing separately) and just focusing on diversification driven by unknown
factors. That is what this function does. It sets up and executes a completely trait-free
version of a HiSSE model.

Thus, all that is required is a tree. The model allows up to 26 possible hidden states in
diversification (denoted by A-Z). Transitions among hidden states are governed by a global
transition rate, q. A "shift" in diversification denotes a lineage tracking some
unobserved, hidden state. An interesting byproduct of this assumption is that distantly
related clades can actually share the same discrete set of diversification parameters.

Note that "hidden state" is a shorthand. We do not mean that there is a single, discrete
character that is solely driving diversification differences. There is some heritable
"thing" that affects rates: but this could be a combination of body size, oxygen
concentration, trophic level, and how many other species are competing in an area.
This is true for HiSSE, but is especially important to grasp for MiSSE. It could
be that there is some single discrete trait that drives everything; it's more likely
that a whole range of factors play a role, and we just slice them up into discrete
categories, the same way we slice up mammals into carnivore / omnivore / herbivore
or plants into woody / herbaceous when the reality is more continuous.

As with \code{hisse}, we employ a modified optimization procedure. In other words, rather
than optimizing birth and death separately, \code{MiSSE} optimizes orthogonal
transformations of these variables: we let tau = birth+death define "net turnover", and
we let eps = death/birth define the \dQuote{extinction fraction}. This reparameterization
alleviates problems associated with overfitting when birth and death are highly
correlated, but both matter in explaining the diversity pattern.

For the \dQuote{root.type} option, we are currently maintaining the previous default of
\dQuote{madfitz}. However, it was recently pointed out by Herrera-Alsina et al. (2018)
that at the root, the individual likelihoods for each possible state should be conditioned
prior to averaging the individual likelihoods across states. This can be set doing
\dQuote{herr_als}. It is unclear to us which is exactly correct, but it does seem that
both \dQuote{madfitz} and \dQuote{herr_als} behave exactly as they should in the case of
character-independent diversification (i.e., reduces to likelihood of tree + likelihood of
trait model). We've also tested the behavior and the likelihood differences are very
subtle and the parameter estimates in simulation are nearly indistinguishable from the
\dQuote{madfitz} conditioning scheme. We provide both options and encourage users to try
both and let us know conditions in which the result vary dramatically under the two root
implementations. We suspect they do not.

Also, note, that in the case of \dQuote{root.type=user} and \dQuote{root.type=equal} are
no longer explicit \dQuote{root.type} options. Instead, either \dQuote{madfitz} or
\dQuote{herr_als} are specified and the \dQuote{root.p} can be set to allow for custom
root options.
}
\value{
\code{MiSSE} returns an object of class \code{misse.fit}. This is a list with
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
\item{$max.tol}{relative optimization tolerance.}
\item{$starting.vals}{The starting values for the optimization.}
\item{$upper.bounds}{the vector of upper limits to the optimization search.}
\item{$lower.bounds}{the vector of lower limits to the optimization
  search.}
\item{$ode.eps}{The ode.eps value used for the estimation.}
\item{$turnover}{The turnover vector used.}
\item{$eps}{The eps vector used.}
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
