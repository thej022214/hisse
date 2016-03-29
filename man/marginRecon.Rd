\name{MarginRecon}
\alias{MarginRecon}
\title{Ancestral State Estimation based on Marginal Reconstruction}
\description{Estimates the likeliest states for both internal nodes and tips of a phylogeny using the marginal reconstruction algorithm.}
\usage{
MarginRecon(phy, data, f, pars, hidden.states=TRUE, four.state.null=FALSE, 
timeslice=NULL, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, 
aic=NULL, verbose=TRUE, n.cores=NULL)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a data matrix containing species information (see Details).}
\item{f}{vector of length 2 with the estimated proportion of extant species in state 0 and 1 that are included in the phylogeny. A value of c(0.25, 0.5) means that 25 percent of species in state 0 and 50 percent of species in state 1 are included in the phylogeny. By default all species are assumed to be known.}
\item{pars}{vector containing the MLE of the parameters.}
\item{hidden.states}{a logical indicating whether the model includes a hidden state. The default is \code{FALSE}.} 
\item{four.state.null}{a logical indicating whether the model is the null-four model. The default is \code{FALSE}.} 
\item{timeslice}{a user-supplied time to split the tree.}
\item{condition.on.survival}{a logical indicating whether the likelihood should be conditioned on the survival of two lineages and the speciation event subtending them (Nee et al. 1994). The default is \code{TRUE}.}
\item{root.type}{indicates whether root prior assumption should based the procedure described by FitzJohn et al. 2009, \dQuote{madfitz}, assumed equal, \dQuote{equal}, or set to user, \dQuote{user}.}
\item{root.p}{a vector indicating fixed root state probabilities. The default is \code{NULL}.}
\item{aic}{the AIC for the model being used for the reconstruction. This is used by the plotting function. The default is \code{NULL}.}
\item{verbose}{a logical indicating whether progress should be printed to screen. The default is \code{TRUE}.}
\item{n.cores}{specifies the number of independent processors to conduct the analysis.. The default is \code{NULL}.}
}
\details{
In this implementation the marginal probability of state i for a focal node is simply the overall likelihood of the tree and data when the state of the focal node is fixed in state i. Note that the likeliest tip states can also be estimated: we observe state 1, but the underlying state could either be state 1A or 1B. Thus, for any given node or tip we traverse the entire tree as many times as there are states in the model. As the size of the tree grows, however, these repeated tree traversals can slow the calculation down considerably. For this reason, we allow the marginal calculation to be conducted in parallel across any number of independent computer processors.

For user-specified \dQuote{root.p}, you should specify the probability for each state. If you are doing a hidden model, there will be four states: 0A, 1A, 0B, 1B. So if you wanted to say the root had to be state 0, you would specify \dQuote{root.p = c(0.5, 0, 0.5, 0)}. 

See help for \dQuote{hisse} for more on other parameters for this function.

}
\value{
\code{MarginRecon} returns an object of class \code{hisse.states}. This is a list with elements:
\item{$node.mat}{the marginal probabilities calculated for each node. They are ordered based on the elements in the edge matrix in dQuote{phylo} format.}
\item{$tip.mat}{the marginal probabilities calculated for each tip. They are ordered based on the order of tip labels in the tree.}
\item{$rate.mat}{a matrix that details the rates for each state combination. This is used by the plotting function.}
\item{$phy}{a phylogenetic tree in the dQuote{phylo} format that contains the states with the highest marginal probability at each internal node.}
}
\examples{
## Not run
# library(diversitree)
# pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
# set.seed(4) 
# phy <- tree.bisse(pars, max.t=30, x0=0)
# sim.dat <- data.frame(names(phy$tip.state), phy$tip.state)
## Now fit HiSSE equivalent with a hidden state for state 1:
# trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
# trans.rates.hisse <- ParDrop(trans.rates.hisse, c(2,3,5,7,8,9,10,12))
# pp.hisse <- hisse(phy, sim.dat, hidden.states=TRUE, turnover.anc=c(1,2,0,3), 
# eps.anc=c(1,2,0,3), trans.rate=trans.rates.hisse)
##Now reconstruct the likeliest states under this model:
# pp <- MarginRecon(phy, sim.dat, f=c(1,1), pars=pp.hisse$solution, 
# hidden.states=TRUE)
}
\references{
Beaulieu, J.M, and B.C. O'Meara. In revision. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. In revision. 

FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Syst. Biol. 58:595-611.

Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary characters' effect on speciation and extinction. Syst. Biol. 56:701-710.

Nee S., May R.M., and Harvey P.H. 1994. The reconstructed evolutionary process. Philos. Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
}
\author{Jeremy M. Beaulieu}
\keyword{models}
