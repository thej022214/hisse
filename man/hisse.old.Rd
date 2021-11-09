\name{hisse.old}
\alias{hisse.old}
\title{Original Hidden State Speciation and Extinction}
\description{Sets up and executes the original HiSSE model (Hidden State Speciation and Extinction) on a phylogeny and character distribution.}
\usage{
hisse.old(phy, data, f=c(1,1), hidden.states=TRUE, turnover.anc=c(1,1,0,0), 
eps.anc=c(1,1,0,0), trans.rate=NULL, turnover.beta=c(0,0,0,0), 
eps.beta=c(0,0,0,0), timeslice=NULL, condition.on.survival=TRUE, 
root.type="madfitz", root.p=NULL, output.type="turnover", sann=TRUE,
sann.its=1000, bounded.search=TRUE, max.tol=.Machine$double.eps^.50,
starting.vals=NULL, turnover.upper=10000, eps.upper=3, trans.upper=100,
ode.eps=0)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a data matrix containing species information (see Details).}
\item{f}{vector of length 2 with the estimated proportion of extant species in state 0 and 1 that are included in the phylogeny. A value of c(0.25, 0.5) means that 25 percent of species in state 0 and 50 percent of species in state 1 are included in the phylogeny. By default all species are assumed to be sampled.}
\item{hidden.states}{a logical indicating whether the model includes a hidden state. The default is \code{FALSE}.} 
\item{turnover.anc}{a vector of length 4, indicating the free parameters associated with the net turnover rates. Default settings is a BiSSE model with fixed turnover rates for both observed states (see Details).} 
\item{eps.anc}{a vector of length 4, indicating the free parameters associated with the extinction fractions. Default settings is a BiSSE model with fixed extinction fractions for both observed states (see Details).} 
\item{trans.rate}{provides the transition rate model.}
\item{turnover.beta}{a vector of length 4, indicating the free parameters associated with time-varying net turnover rates (see Details).}
\item{eps.beta}{a vector of length 4, indicating the free parameters associated with time-varying extinction fractions (see Details).}
\item{timeslice}{a user-supplied time to split the tree.}
\item{condition.on.survival}{a logical indicating whether the likelihood should be conditioned on the survival of two lineages and the speciation event subtending them (Nee et al. 1994). The default is \code{TRUE}.}
\item{root.type}{indicates whether root summarization follow the procedure described by FitzJohn et al. 2009, \dQuote{madfitz} or Herrera-Alsina et al. 2018, \dQuote{herr_als}.}
\item{root.p}{a vector indicating fixed root state probabilities. The default is \code{NULL}.}
\item{output.type}{indicates whether the rates should be printed onscreen as the optimized variables, \dQuote{turnover}, transformed to reflect net diversification, \dQuote{net.div}, or transformed to reflect \eqn{\lambda}{lambda} and \eqn{\mu}{mu}, \dQuote{raw}.}
\item{sann}{a logical indicating whether a two-step optimization procedure is to be used. The first includes a simulate annealing approach, with the second involving a refinement using \code{subplex}. The default is \code{TRUE}.}
\item{sann.its}{a numeric indicating the number of times the simulated annealing algorithm should call the objective function.}
\item{bounded.search}{a logical indicating whether or not bounds should be enforced during optimization. The default is is \code{TRUE}.}
\item{max.tol}{supplies the relative optimization tolerance to \code{subplex}.}
\item{starting.vals}{a vector of starting values to be used instead of the default settings. These are just three values given in the following order: turnover (1), extinction fraction (2), and a single transition rate (3)}
\item{turnover.upper}{sets the upper bound for the turnover parameters. The default upper bound assumes an event occurs every 100 years.}
\item{eps.upper}{sets the upper bound for the extinction fraction parameters.}
\item{trans.upper}{sets the upper bound for the transition rate parameters.}
\item{ode.eps}{sets the tolerance for the integration at the end of a branch. Essentially if the sum of compD is less than this tolerance, then it assumes the results are unstable and discards them. The default is set to zero, but in testing a value of 1e-8 can sometimes produce stable solutions for both easy and very difficult optimization problems.}
}
\details{
This function sets up and executes the original HiSSE model. The model closely follows \code{diversitree}, although here we employ modified optimization procedures. For example, rather than optimizing birth and death separately, hisse optimizes orthogonal transformations of these variables: we let tau = birth+death define "net turnover", and we let eps = death/birth define the \dQuote{extinction fraction}. This reparameterization alleviates problems associated with overfitting when birth and death are highly correlated, but both matter in explaining the diversity pattern. As for data file format, \code{hisse} expects a two column matrix or data frame, with the first column containing the species names and the second containing the binary character information. Note that the order of the data file and the names in the \dQuote{phylo} object need not be in the same order; \code{hisse} deals with this internally. Also, the character information should be binary and coded as 0 and 1, otherwise the function will misbehave. However, if the state for a species is unknown, a user can specify this with a 2, and the state will be considered maximally ambiguous.

To setup a model, users input vectors containing values to indicate how many free parameters are to be estimated for each of the variables in the model. For example, the \dQuote{turnover.anc} input vector is set by default as c(1,1,0,0). This means for state 0 and state 1, we are allowing one free parameter to define the net turnover rate (birth+death) in the model. This is essentially a BiSSE model with fixed turnover rates. Now, say we want to include separate turnover rates for both states we would simply input c(1,2,0,0). The last two entries, which in the preceding example are set to zero, correspond to the hidden states; the third entry corresponds to a hidden state associated with observed state 0, such that 0A (hidden state absent) is the first entry, and 0B (hidden state present) is the third entry. So, to set up a model with three turnover rates, where we include a free parameter for a hidden state associated with state 0 we input c(1,2,3,0). A full model would thus be c(1,2,3,4), which corresponds to four separate net turnover rates, for states 0A (first entry), 1A (second entry), 0B (third entry), and 1B (fourth entry). Extinction fraction, or \dQuote{eps.anc}, follows the same format, though including a zero for a state we want to include in the model corresponds to no extinction, which is the Yule equivalent. In general, we follow this format to make it easier to generate a large set of nested models. Once the model is specified, the parameters can be estimated using the subplex routine (default), or use a two-step process (i.e., sann=TRUE) that first employs a stochastic simulated annealing procedure, which is later refined using the subplex routine.

The \dQuote{trans.rate} input is the transition model and has an entirely different setup than turnover and extinction rates. See \code{TransMatMaker} function for more details. 

For user-specified \dQuote{root.p}, you should specify the probability for each state. If you are doing a hidden model, there will be four states: 0A, 1A, 0B, 1B. So if you wanted to say the root had to be state 0, you would specify \dQuote{root.p = c(0.5, 0, 0.5, 0)}.

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

Finally, the options \dQuote{.beta} and \dQuote{timeslice} are included, but neither have been tested -- needless to say, use at your own risk (but really, though, you should probably forget that these options exist for the time being). The \dQuote{.beta} provides a means for testing for time-varying rates, whereas \dQuote{timeslice} splits the tree to allow the process to vary before and after some user defined time period. These options will be further developed in due course.
}
\value{
\code{hisse.old} returns an object of class \code{hisse.fit}. This is a list with elements:
\item{$loglik}{the maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample-size.}
\item{$solution}{a matrix containing the maximum likelihood estimates of the model parameters.}
\item{$index.par}{an index matrix of the parameters being estimated.}
\item{$f}{user-supplied sampling frequencies.}
\item{$hidden.states}{a logical indicating whether a hidden state was included in the model.}
\item{$condition.on.surivival}{a logical indicating whether the likelihood was conditioned on the survival of two lineages and the speciation event subtending them.}
\item{$root.type}{indicates the user-specified root prior assumption.}
\item{$root.p}{indicates whether the user-specified fixed root probabilities.}
\item{$timeslice}{indicates whether the user-specified timeslice that split the tree.}
\item{$phy}{user-supplied tree}
\item{$data}{user-supplied dataset}
\item{$iterations}{number of iterations of the likelihood search that were executed.}
\item{$output.type}{the user-specified output.type to be printed on the screen.}
\item{$max.tol}{relative optimization tolerance.}
\item{$upper.bounds}{the vector of upper limits to the optimization search.}
\item{$lower.bounds}{the vector of lower limits to the optimization search.}
}
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 

FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Syst. Biol. 58:595-611.

Herrera-Alsina, L., P. van Els, and R.S. Etienne. 2018. Detecting the dependence of diversification on multiples traits from phylogenetic trees and trait data. Systematic Biology, 68:317-328.

Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary characters effect on speciation and extinction. Syst. Biol. 56:701-710.

Nee S., May R.M., and Harvey P.H. 1994. The reconstructed evolutionary process. Philos. Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
}

\author{Jeremy M. Beaulieu}
\keyword{models}
