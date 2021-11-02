\name{TipCorrelation}
\alias{TipCorrelation}
\title{Phylogenetic independent contrasts using tip rates}
\description{Performs linear regression between phylogenetic independent contrasts (PICs) of tip rates and continuous trais}
\usage{
TipCorrelation(phy, tip.rate, trait, log=TRUE, remove.cherries=TRUE, 
scaled=TRUE, positivise=TRUE, use.lmorigin=TRUE)
}       
\arguments{
\item{phy}{an ultrametric phylogenetic tree of class 'phylo'.} 
\item{tip.rate}{an object of class 'named numeric' containing tip.rates}
\item{trait}{an object of class 'named numeric' containing continuous trait data}
\item{log}{Should tip.rate and trait values be logged before regression? The default is TRUE.}
\item{remove.cherries}{Should PICs of 'cherries' be pruned before regression? The default is TRUE.}
\item{scaled}{Should PICs be scaled with their expected variances? The default to TRUE.}
\item{positivise}{Should PICs be positivised before regression? The default is TRUE.}
\item{use.lmorigin}{Should a regression-through-origin be performed instead of regular linear model? The default is TRUE.}
}
\details{
Tip rates as those obtained with \code{MiSSEGreedy} are analogous to a continuous trait evolving in the tree and must also be corrected for phylogenetic 
non-independence before regression analyses. This function does that by using PICs (Felsenstein, 1985) and gives the additional option of removing PICs 
from 'cherries' from analyses. Cherries (theoretically) inherit the exact same rate class probabilities in any model that uses just branch lengths to estimate 
tip rates. For that reason, they may: (1) present identical tip rates, forcing the slope of the regression to be close to 0 since all PICs for cherries will be 0; 
and (2) constitute pseudoreplicates in the analyses. We suspect that, for that reason, it may make sense to prune them out from any PIC analyses that uses tip-rate
correlations from tree-only diversification approaches (see Vasconcelos et al. in prep.).
}
\value{
Returns an object that contains:
\item{$correlation}{summary of correlation results.}
\item{$tip.rate PIC}{tip rate PICs for each node.}
\item{$trait PIC}{trait PICs for each node.}
}
\references{
Vasconcelos, T, B.C. O'Meara, and J.M. Beaulieu. In prep.
Felsenstein, J. (1985). Phylogenies and the comparative method. The American Naturalist, 125, 1-15.
}
\author{Thais Vasconcelos and Brian O'Meara}
