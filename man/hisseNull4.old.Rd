\name{hisse.null4.old}
\alias{hisse.null4.old}
\title{Original Four state trait-independent Hidden State Speciation and Extinction}
\description{Sets up and executes the original four state trait-independent HiSSE model (Hidden State Speciation and Extinction) on a phylogeny and character set.}
\usage{
hisse.null4.old(phy, data, f=c(1,1), turnover.anc=rep(c(1,2,3,4),2), 
eps.anc=rep(c(1,2,3,4),2), trans.type="equal", condition.on.survival=TRUE, 
root.type="madfitz", root.p=NULL,  output.type="turnover", sann=TRUE, 
sann.its=1000, bounded.search=TRUE, max.tol=.Machine$double.eps^.50,
starting.vals=NULL, turnover.upper=10000, eps.upper=3, trans.upper=100,
ode.eps=0)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a data matrix containing species information (see Details).}
\item{f}{vector of length 2 with the estimated proportion of extant species in state 0 and 1 that are included in the phylogeny. A value of c(0.25, 0.5) means that 25 percent of species in state 0 and 50 percent of species in state 1 are included in the phylogeny. By default all species are assumed to be sampled.}
\item{turnover.anc}{a vector of length 8, indicating the free parameters associated with the net turnover rates. Default setting assumes character independent diversification (see Details).} 
\item{eps.anc}{a vector of length 8, indicating the free parameters associated with the extinction fractions. Default setting assumes character independent diversification (see Details).} 
\item{trans.type}{provides the type of transition rate model. Currently this model allows two types: \dQuote{equal}, the default, which assumes all transitions are equal, and \dQuote{three.rate}, that assumes three rates (see Details).} 
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
This function sets up and executes a four-state trait independent HiSSE model. The model closely follows \code{hisse.old}. However, note that this function is no longer necessary and can be called and evaluated directly using the new \code{hisse} function.

Like \code{hisse.old}, users input vectors containing values to indicate how many free parameters are to be estimated for each of the variables in the model. However, the null-four model assumes that \dQuote{turnover.anc} and \dQuote{eps.anc} are linked between the two observed states. Thus, users are unlikely to alter the inputs much, aside from perhaps fixing \dQuote{turnover.anc} or \dQuote{eps.anc} to be equal across the four hidden states, where the \dQuote{turnover.anc} input vector is set as rep(c(1,1,1,1),2). For a Yule equivalent, the input vector for \dQuote{eps.anc} would be rep(c(0,0,0,0),2). For how to setup a null-two model see the example code below.

For user-specified \dQuote{root.p}, you should specify the probability for each state. See help for \dQuote{hisse.old} for more on other parameters for this function.

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
\code{hisse.null4.old} returns an object of class \code{hisse.fit}. This is a list with elements:
\item{$loglik}{the maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample-size.}
\item{$solution}{a matrix containing the maximum likelihood estimates of the model parameters.}
\item{$index.par}{an index matrix of the parameters being estimated.}
\item{$f}{user-supplied sampling frequencies.}
\item{$condition.on.surivival}{a logical indicating whether the likelihood was conditioned on the survival of two lineages and the speciation event subtending them.}
\item{$root.type}{indicates the user-specified root prior assumption.}
\item{$root.p}{indicates whether the user-specified fixed root probabilities.}
\item{$phy}{user-supplied tree}
\item{$data}{user-supplied dataset}
\item{$output.type}{the user-specified output.type to be printed on the screen.}
\item{$trans.type}{the user-specified transition model.}
\item{$trans.mat}{the index matrix that specifies the free parameters in the transition model.} 
\item{$max.tol}{relative optimization tolerance.}
\item{$upper.bounds}{the vector of upper limits to the optimization search.}
\item{$lower.bounds}{the vector of lower limits to the optimization search.}
}
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 
}

\author{Jeremy M. Beaulieu}
\keyword{models}
