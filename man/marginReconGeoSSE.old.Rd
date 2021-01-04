\name{MarginReconGeoSSE.old}
\alias{MarginReconGeoSSE.old}
\title{Original Ancestral State Estimation based on Marginal Reconstruction for
  the GeoSSE and GeoHiSSE models.}
\description{Estimates the likeliest states for both internal nodes and tips of a phylogeny using the marginal reconstruction algorithm.}
\usage{
MarginReconGeoSSE.old(phy, data, f, pars, hidden.areas=TRUE, 
assume.cladogenetic=TRUE, condition.on.survival=TRUE, 
root.type="madfitz", root.p=NULL, AIC=NULL, verbose=TRUE, 
n.cores=NULL)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a data matrix containing species information (see Details).}
\item{f}{vector of length 2 with the estimated proportion of extant species in state 0 and 1 that are included in the phylogeny. A value of c(0.25, 0.5) means that 25 percent of species in state 0 and 50 percent of species in state 1 are included in the phylogeny. By default all species are assumed to be known.}
\item{pars}{vector containing the MLE of the parameters.}
\item{hidden.areas}{a logical indicating whether the model includes hidden areas. The default is \code{FALSE}.} 
\item{assume.cladogenetic}{a logical indicating whether the model
  assumes that changes occurr at the nodes. The default is
  \code{TRUE}. Setting this to \code{FALSE} will result in a different
  models from GeoSSE. Please check \code{hisse::GeoHiSSE.old} for more details.}
\item{condition.on.survival}{a logical indicating whether the likelihood should be conditioned on the survival of two lineages and the speciation event subtending them (Nee et al. 1994). The default is \code{TRUE}.}
\item{root.type}{indicates whether root prior assumption should based the procedure described by FitzJohn et al. 2009, \dQuote{madfitz}, assumed equal, \dQuote{equal}, or set to user, \dQuote{user}.}
\item{root.p}{a vector indicating fixed root state probabilities. The default is \code{NULL}.}
\item{AIC}{the AIC for the model being used for the reconstruction. This is used by the plotting function. The default is \code{NULL}.}
\item{verbose}{a logical indicating whether progress should be printed to screen. The default is \code{TRUE}.}
\item{n.cores}{specifies the number of independent processors to conduct the analysis.. The default is \code{NULL}.}
}
\details{
In this implementation the marginal probability of area i for a focal node is simply the overall likelihood of the tree and data when the area of the focal node is fixed in area i. Note that the likeliest tip areas can also be estimated: we observe area 1, but the underlying area could either be area 1A or 1B. Thus, for any given node or tip we traverse the entire tree as many times as there are areas in the model. As the size of the tree grows, however, these repeated tree traversals can slow the calculation down considerably. For this reason, we allow the marginal calculation to be conducted in parallel across any number of independent computer processors.

For user-specified \dQuote{root.p}, you should specify the probability
for each area. If you are doing a hidden model, the number of parameters
will depend on the number of hidden areas included. For a two classes
model there are six areas: A0, B0, AB0, A1, B1, AB1. So if you wanted to
say the root had to be area A, you would specify \dQuote{root.p = c(0.5,
  0, 0, 0.5, 0, 0)}. There is 'root area is A and there is an equal
chance for hidden area A0 or A1'.

See help for \dQuote{GeoHiSSE.old} for more on other parameters for this function.

}
\value{
\code{MarginReconGeoSSE.old} returns an object of class \code{geohisse.states}. This is a list with elements:
\item{$node.mat}{the marginal probabilities calculated for each node. They are ordered based on the elements in the edge matrix in \dQuote{phylo} format.}
\item{$tip.mat}{the marginal probabilities calculated for each tip. They are ordered based on the order of tip labels in the tree.}
\item{$rate.mat}{a matrix that details the rates for each state combination. This is used by the plotting function.}
\item{$phy}{a phylogenetic tree in the \dQuote{phylo} format that contains the states with the highest marginal probability at each internal node.}
}
\examples{
\donttest{
## Need some examples here.
}
}
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 

FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Syst. Biol. 58:595-611.

Goldberg, E. E., L. T. Lancaster, and R. H. Ree. 2011. Phylogenetic Inference of Reciprocal Effects between Geographic Range Evolution and Diversification. Syst. Biol. 60:451-465.

Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary characters' effect on speciation and extinction. Syst. Biol. 56:701-710.

Nee S., May R.M., and Harvey P.H. 1994. The reconstructed evolutionary process. Philos. Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
}
\author{Jeremy M. Beaulieu}
\keyword{models}
