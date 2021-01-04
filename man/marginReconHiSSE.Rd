\name{MarginReconHiSSE}
\alias{MarginReconHiSSE}
\title{Ancestral State Estimation based on Marginal Reconstruction for
  the HiSSE model.}
\description{Estimates the likeliest states for both internal nodes and tips of a phylogeny using the marginal reconstruction algorithm.}
\usage{
MarginReconHiSSE(phy, data, f, pars, hidden.states=1, 
condition.on.survival=TRUE,  root.type="madfitz", root.p=NULL, 
AIC=NULL, get.tips.only=FALSE, verbose=TRUE, n.cores=NULL, dt.threads=1)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a matrix (or dataframe) with three columns. The first column 
  containing the species names and the second and third containing the binary 
  character information. Character "0" is on column 2 and chracter "1" is on 
  column 3. A value of 0 means character absent and a value of 1 character 
  present. The input of data follows a Pagel model. See 'Details'.}
\item{f}{vector of length 2 with the estimated proportion of extant species in 
  state 0 and 1 that are included in the phylogeny. A value of c(0.25, 0.5) means 
  that 25 percent of species in state 0 and 50 percent of species in state 1 are 
  included in the phylogeny. By default all species are assumed to be sampled.}
\item{pars}{vector containing the MLE of the parameters.}
\item{hidden.states}{a numeric indicating the number of shifts. The default is \code{1}.} 
\item{condition.on.survival}{a logical indicating whether the likelihood should be conditioned on the survival of two lineages and the speciation event subtending them (Nee et al. 1994). The default is \code{TRUE}.}
\item{root.type}{indicates whether root prior assumption should based the procedure described by FitzJohn et al. 2009, \dQuote{madfitz}, assumed equal, \dQuote{equal}, or set to user, \dQuote{user}.}
\item{root.p}{a vector indicating fixed root state probabilities. The default is \code{NULL}.}
\item{AIC}{the AIC for the model being used for the reconstruction. This is used by the plotting function. The default is \code{NULL}.}
\item{get.tips.only}{a logical indicating whether just tip reconstructions should be output. The default is \code{FALSE}.}
\item{verbose}{a logical indicating whether progress should be printed to screen. The default is \code{TRUE}.}
\item{n.cores}{specifies the number of independent processors to conduct the analysis.. The default is \code{NULL}.}
\item{dt.threads}{sets the number of threads available to data.table. In practice this need not change from the default of 1 thread, as we have not seen any speedup from allowing more threads.}
}
\details{
This is the marginal reconstruction algorithm for the newer, faster version of HiSSE. In this implementation the marginal probability of state i for a focal node is simply the overall likelihood of the tree and data when the state of the focal node is fixed in state i. Note that the likeliest tip states can also be estimated: we observe state 1, but the underlying state could either be state 1A or 1B. Thus, for any given node or tip we traverse the entire tree as many times as there are states in the model. As the size of the tree grows, however, these repeated tree traversals can slow the calculation down considerably. For this reason, we allow the marginal calculation to be conducted in parallel across any number of independent computer processors.

See help for \dQuote{hisse} for more on other parameters for this function.

}
\value{
\code{MarginReconHiSSE} returns an object of class \code{hisse.states}. This is a list with elements:
\item{$node.mat}{the marginal probabilities calculated for each node. They are ordered based on the elements in the edge matrix in \dQuote{phylo} format.}
\item{$tip.mat}{the marginal probabilities calculated for each tip. They are ordered based on the order of tip labels in the tree.}
\item{$rate.mat}{a matrix that details the rates for each state combination. This is used by the plotting function.}
\item{$phy}{a phylogenetic tree in the \dQuote{phylo} format that contains the states with the highest marginal probability at each internal node.}
}
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 

FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Syst. Biol. 58:595-611.

Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary characters' effect on speciation and extinction. Syst. Biol. 56:701-710.

Nee S., May R.M., and Harvey P.H. 1994. The reconstructed evolutionary process. Philos. Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
}
\author{Jeremy M. Beaulieu}
\keyword{models}
