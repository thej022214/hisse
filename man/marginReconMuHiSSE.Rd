\name{MarginReconMuHiSSE}
\alias{MarginReconMuHiSSE}
\title{Ancestral State Estimation based on Marginal Reconstruction for
  the MuSSE and MuHiSSE models.}
\description{Estimates the likeliest states for both internal nodes and tips of a phylogeny using the marginal reconstruction algorithm.}
\usage{
MarginReconMuHiSSE(phy, data, f, pars, hidden.states=TRUE, 
condition.on.survival=TRUE,  root.type="madfitz", root.p=NULL, 
aic=NULL, verbose=TRUE, n.cores=NULL)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a matrix (or dataframe) with three columns. The first column 
  containing the species names and the second and third containing the binary 
  character information. Character "0" is on column 2 and chracter "1" is on 
  column 3. A value of 0 means character absent and a value of 1 character 
  present. The input of data follows a Pagel model. See 'Details'.}
\item{f}{vector of length 4 with the estimated proportion of extant
  species in 00, 01, 10, and 11 that are included in the phylogeny. A value of 
  c(0.50, 0.25, 0.125, 0.125) means that 50 percent of species in combination '00', 25 
  percent in '01' and 12.5 percent in '10' and '11'. By default all species are assumed 
  to be sampled.}
\item{pars}{vector containing the MLE of the parameters.}
\item{hidden.states}{a logical indicating whether the model includes hidden states. The default is \code{FALSE}.} 
\item{condition.on.survival}{a logical indicating whether the likelihood should be conditioned on the survival of two lineages and the speciation event subtending them (Nee et al. 1994). The default is \code{TRUE}.}
\item{root.type}{indicates whether root prior assumption should based the procedure described by FitzJohn et al. 2009, \dQuote{madfitz}, assumed equal, \dQuote{equal}, or set to user, \dQuote{user}.}
\item{root.p}{a vector indicating fixed root state probabilities. The default is \code{NULL}.}
\item{aic}{the AIC for the model being used for the reconstruction. This is used by the plotting function. The default is \code{NULL}.}
\item{verbose}{a logical indicating whether progress should be printed to screen. The default is \code{TRUE}.}
\item{n.cores}{specifies the number of independent processors to conduct the analysis.. The default is \code{NULL}.}
}
\details{
In this implementation the marginal probability of state combination i for a focal node is simply the overall likelihood of the tree and data when the state of the focal node is fixed in state combination i. Note that the likeliest tip states can also be estimated: we observe state 00, but the underlying state could either be 00A or 00B. Thus, for any given node or tip we traverse the entire tree as many times as there are states in the model. As the size of the tree grows, however, these repeated tree traversals can slow the calculation down considerably. For this reason, we allow the marginal calculation to be conducted in parallel across any number of independent computer processors.

For user-specified \dQuote{root.p}, you should specify the probability
for each state combination. If you are doing a hidden model, the number of parameters
will depend on the number of hidden states included. For a two classes
model there are eight states: 00A, 01A, 10A, 11A, 00B, 01B, 10B, and 11B. So if you wanted to
say the root had to be in state 00, you would specify \dQuote{root.p = c(0.5,
  0, 0, 0, 0.5, 0, 0, 0, 0)}. There is 50 percent chance the root state is 00 and there is an equal
chance for hidden state A or B.

See help for \dQuote{MuHiSSE} for more on other parameters for this function.

}
\value{
\code{MarginRecon} returns an object of class \code{muhisse.states}. This is a list with elements:
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

Nakov, T., Beaulieu, J.M., and Alverson, A.J. 2018. Freshwater diatoms diversify faster than marine in both planktonic and benthic habitats. bioRxiv, doi: https://doi.org/10.1101/406165.

Nee S., May R.M., and Harvey P.H. 1994. The reconstructed evolutionary process. Philos. Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
}
\author{Jeremy M. Beaulieu}
\keyword{models}
