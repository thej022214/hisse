\name{MiSSEGreedy}
\alias{MiSSEGreedy}
\title{Character-free State Speciation and Extinction searching greedily}
\description{Sets up and executes a set of MiSSE models (Missing State Speciation
 and Extinction) on a phylogeny, varying the number of parameters for turnover and extinction fraction and stopping when models stop being very good.}
\usage{
MiSSEGreedy(phy, f=1, turnover.tries=sequence(26), eps.constant=c(TRUE,FALSE), 
stop.count=2, stop.deltaAICc=10, condition.on.survival=TRUE, root.type="madfitz",
root.p=NULL, sann=FALSE, sann.its=10000, bounded.search=TRUE, 
max.tol=.Machine$double.eps^.50, starting.vals=NULL, turnover.upper=10000, eps.upper=3, 
trans.upper=100, restart.obj=NULL, ode.eps=0)
}
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and
  with internal nodes labeled denoting the ancestral selective regimes.}
\item{f}{the estimated proportion of extant species included in the phylogeny. A value of
  0.50 means that 50 percent of species are contained in the. By default all species are assumed
  to be sampled.}
\item{turnover.tries}{a numeric vector of how many turnover categories to try. See 'Details'.}
\item{eps.constant}{a boolean vector of whether to have eps be constant (one parameter if TRUE) or equal in complexity to the number of turnover categories. See 'Details'.}
\item{stop.count}{how many steps to take of models that have a delta AICc worse than stop.deltaAICc}
\item{stop.deltaAICc}{how bad compared to the best does a model have to be to far enough outside we stop looking at even more complex ones?}
\item{condition.on.survival}{a logical indicating whether the likelihood
  should be conditioned on the survival of two lineages and the
  speciation event subtending them (Nee et al. 1994). The default is \code{TRUE}.}
\item{root.type}{indicates whether root summarization follow the procedure described by FitzJohn et al. 2009, \dQuote{madfitz} or Herrera-Alsina et al. 2018, \dQuote{herr_als}.}
\item{root.p}{a vector indicating fixed root state probabilities. The
  default is \code{NULL}.}
\item{sann}{a logical indicating whether a two-step optimization
  procedure is to be used. The first includes a simulate annealing
  approach, with the second involving a refinement using
  \code{subplex}. The default is \code{FALSE}.}
\item{sann.its}{a numeric indicating the number of times the simulated
  annealing algorithm should call the objective function.}
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
}
\details{
See the \code{MiSSE} function for description of the method overall. It requires a set number of hidden state categories, but finding the best number of categories can be hard. For example, one could 1 to 26 different turnover categories and 1 to 26 possible extinction fraction categories. For most cases, we suspect that it makes sense to have the number of extinction fraction categories either equal to the number of turnover categories or set to the same category over the tree, but there are actually a lot of possibilities: have turnover=c(1,2,3) and eps=c(1,2,1), for example, or turnover=c(1,1,1) and eps=c(1,2,3). As a shortcut, this function tries working up in complexity. By default, it starts with one turnover rate and one extinction fraction. Then, it tries two of each, then three of each, and so forth. Then it starts over with increasing number of turnover rates but just keeping one extinction fraction. This is 52 possible models, which, even with MiSSE's speed [that's a joke] can take a fair bit of time. It's thus possible to set a stopping criterion so that once models are complex enough they're not going to have a lot of AICc weight, the analysis stops. This is NOT where the models stop being significant (if you're looking for significance, note you don't get that from AICc), but where they probably will not contribute much to the model averaged parameter estimates and so may not be worth the bother of searching further. However, there's no guarantee that this is a wise decision: you could stop with 10 turnover rates, and 11 and 12 are far worse for AICc, but it could be that 13 turnover rates are much better than anything else (for that matter, the best number of turnover rates could be 42, even though MiSSe's current code can only go up to 26). It could be that turnover only needs five hidden states but extinction fraction works best with nine. And of course, in reality, the truth is that there is an infinite set of hidden rate parameters: a passing cloud blocks a tiny bit of energy for photosynthesis, so for that moment the rate of extinction for plants underneath the shadow is a tiny bit higher. You should not be using this to test hypotheses about how many hidden factors there are, but rather as a way to get a good enough model to estimate rates on a tree. Note that if allowed to use both eps settings (same number of parameters as turnover or one parameter) it concatenates all the models it examines but it stops separately for each eps setting, allowing the number of free parameters to increase smoothly.

By default, this will try from 1 to 26 different turnover states, with either 1 or n different extinction fractions. If you just want to try up to five states, then turnover.tries=c(1:5). If you want to try 1:5, 10, 15, 20, 25, then turnover.tries=c(1:5, 10, 15, 20, 25). You can change the order in which models are tried, too: turnover.tries=c(16, 14, 12, 20) is possible. If you always want the complexity of the extinction fraction to match that of turnover, eps.constant=FALSE; if you want to have it always one parameter, eps.constant=TRUE; if you want to try both, eps.constant=c(TRUE, FALSE).

You can change how quickly the function stops trying new models with stop.count and stop.deltaAICc. Once it gets a model that is at least stop.deltaAICc worse than the *current* best model, it starts counting; once it gets stop.count of these in a row, it stops. Since this is based on current best AICc, and we start with simple models, there's an asymmetry: a terrible model with no rate variation is always included, but a slightly less terrible model with 26 turnover rates might never be run.
}
\value{
\code{MiSSE} returns a list of class \code{misse.fit} objects.
}
\examples{
\donttest{
  library(ape)
  phy <- rcoal(30)
  result <- MiSSEGreedy(phy)
  print(length(result))
}
}
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models
of trait-dependent speciation and extinction. Syst. Biol. 65:583-601.

FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating trait-dependent speciation
and extinction rates from incompletely resolved phylogenies. Syst. Biol. 58:595-611.

Herrera-Alsina L., van Els P., and Etienne R.S. 2018. Detecting the dependence of
diversification on multiples traits from phylogenetic trees and trait data. Systematic
Biology In press.

Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary characters effect on
speciation and extinction. Syst. Biol. 56:701-710.

Nee S., May R.M., and Harvey P.H. 1994. The reconstructed evolutionary process. Philos.
Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
}

\author{Brian C. O'Meara}
\keyword{models}
