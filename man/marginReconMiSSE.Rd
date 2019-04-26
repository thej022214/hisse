\name{MarginReconMiSSE}
\alias{MarginReconMiSSE}
\title{Ancestral State Estimation based on Marginal Reconstruction for
  the MiSSE models.}
\description{Estimates the likeliest states for both internal nodes and tips of a phylogeny using the marginal reconstruction algorithm.}
\usage{
MarginReconMiSSE(phy, f, pars, hidden.states=2, 
condition.on.survival=TRUE,  root.type="madfitz", root.p=NULL, 
aic=NULL, get.tips.only=FALSE, verbose=TRUE, n.cores=NULL)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes.}
\item{f}{the estimated proportion of extant species included in the phylogeny. A value of 
  0.50 means that 50 percent of species are contained in the. By default all species are assumed 
  to be sampled.}
\item{pars}{vector containing the MLE of the parameters.}
\item{hidden.states}{a numeric indicating the number of shifts. The default is \code{2}.} 
\item{condition.on.survival}{a logical indicating whether the likelihood should be conditioned on the survival of two lineages and the speciation event subtending them (Nee et al. 1994). The default is \code{TRUE}.}
\item{root.type}{indicates whether root prior assumption should based the procedure described by FitzJohn et al. 2009, \dQuote{madfitz}, assumed equal, \dQuote{equal}, or set to user, \dQuote{user}.}
\item{root.p}{a vector indicating fixed root state probabilities. The default is \code{NULL}.}
\item{aic}{the AIC for the model being used for the reconstruction. This is used by the plotting function. The default is \code{NULL}.}
\item{verbose}{a logical indicating whether progress should be printed to screen. The default is \code{TRUE}.}
\item{get.tips.only}{a logical indicating whether just tip reconstructions should be output. The default is \code{FALSE}.}
\item{n.cores}{specifies the number of independent processors to conduct the analysis. The default is \code{NULL}.}
}
\details{
In this implementation the marginal probability of state combination i for a focal node is simply the overall likelihood of the tree and data when the state of the focal node is fixed in state combination i.

See help for \dQuote{MiSSE} for more on other parameters for this function.

}
\value{
\code{MarginReconMiSSE} returns an object of class \code{misse.states}. This is a list with elements:
\item{$node.mat}{the marginal probabilities calculated for each node. They are ordered based on the elements in the edge matrix in \dQuote{phylo} format.}
\item{$tip.mat}{the marginal probabilities calculated for each tip. They are ordered based on the order of tip labels in the tree.}
\item{$rate.mat}{a matrix that details the rates for each state combination. This is used by the plotting function.}
\item{$phy}{a phylogenetic tree in the \dQuote{phylo} format that contains the states with the highest marginal probability at each internal node.}
}
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 

FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating
trait-dependent speciation and extinction rates from incompletely
resolved phylogenies. Syst. Biol. 58:595-611.

Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary characters' effect on speciation and extinction. Syst. Biol. 56:701-710.

Nee S., May R.M., and Harvey P.H. 1994. The reconstructed evolutionary process. Philos. Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
}
\author{Jeremy M. Beaulieu}
\keyword{models}
