\name{generateMiSSEGreedyCombinations}
\alias{generateMiSSEGreedyCombinations}
\title{Combinations for Character-free State Speciation and Extinction Searching Greedily}
\description{Sets up a set of MiSSE models (Missing State Speciation
 and Extinction) on a phylogeny, varying the number of parameters for turnover and extinction fraction.}
\usage{
generateMiSSEGreedyCombinations(max.param=52, turnover.tries=sequence(26), 
eps.tries=sequence(26), fixed.eps.tries=NA, vary.both=TRUE, shuffle.start=TRUE)
}
\arguments{
\item{max.param}{how many parameters to try estimating at the most.}
\item{turnover.tries}{vector of number of free turnover parameters.}
\item{eps.tries}{vector of number of free eps parameters.}
\item{fixed.eps.tries}{fixed eps values to use. By default this is set to NA to allow estimating. However, a vector of values can supplied as fixed value for eps (and an NA to allow estimating, as well).}
\item{vary.both}{if TRUE, allows models that have multiple parameters for turnover and eps; if FALSE, if turnover has multiple hidden states eps can only have one estimated value and vice versa.}
\item{shuffle.start}{if TRUE, instead of generating models strictly by increasing complexity shuffles them around so that simpler models tend to be earlier but with some more complex models mixed in}
}
\details{
This creates the set of combinations of models to run through MiSSEGreedy as a data.frame. It has columns for the number of rates to estimate for turnover, the number of values to estimate for eps, and any fixed values for eps to apply to the whole tree. You can add your own rows or delete some to this data.frame to add more or fewer combinations. By default, this comes out ordered so that simpler models are run first in \code{MiSSEGreedy} but that is not required (but wise for most use cases), and you can reorder if you wish.

It can be worth considering how much information your tree has for diversification rates. A fully resolved tree with N taxa has 2N-2 branches and only N-1 internal node heights, which is more relevant for diversification models.  So on a 50 taxon tree, there are 49 node heights -- it's a bit ambitious to think you can extract, say, 5 different parameters (3 turnover rates, 1 eps value, 1 transition rate) from such a tree. People clearly try to extract more -- some methods claim to use information in a tree to extract a different speciation or net diversification rate for every single tip -- but that could be a tad ambitious. So setting \code{max.param} to a low value makes a lot of sense. \code{max.param = floor(ape::Ntip(phy)/10)} is ridiculously optimistic (dividing by 100 or more is probably more conservative) but if things are running well \code{MiSSEGreedy} should stop before getting to the models that are too complex.

\code{turnover.tries} and \code{eps.tries} set how many turnover and eps hidden states to try, respectively. If you wanted to try only 1, 3, and 7 hidden states for turnover you would set \code{turnover.tries = c(1, 3, 7)} for example.

Estimating extinction rates is hard. This affects all diversification models (even if all you want and look at is speciation rate, extinction rate estimates still affect what this is as they affect the likelihood). It is most noticeable in MiSSE with eps, the extinction fraction (extinction rate divided by speciation rate). One option, following Magallon & Sanderson (2001), is to set extinction fraction at set values. By default, we use theirs, 0 (meaning a Yule model - no extinction) or 0.9 (a lot of extinction, though still less than paleontoligists find). You can set your own in \code{fixed.eps.tries}. If you only want to use fixed values, and not estimate, get rid of the NA, as well. However, don't \dQuote{cheat} -- if you use a range of values for fixed.eps, it's basically doing a search for this, though the default AICc calculation doesn't dQuote{know} this to penalize it for another parameter.

HiSSE and thus MiSSE assume that a taxon has a particular hidden state (though they recognize that there can be uncertainty in which state it actually has). Thus, they're written to assume that we dQuote{paint} these states on the tree and a given state affects both turnover and eps. So if turnover has four hidden states, eps has four hidden states. They can be constrained: the easiest way is to have, say, turnover having an independent rate for each hidden state and eps having the same rate for all the hidden states. If \code{vary.both} is set to FALSE, all models are of this sort: if turnover varies, eps is constant across all hidden states, or vice versa. Jeremy Beaulieu prefers this. If \code{vary.both} is set to TRUE, both can vary: for example, there could be five hidden states for both turnover and eps, but turnover lets each of these have a different rate, but eps only allows three values (so that eps_A and eps_D might be forced to be equal, and eps_B and eps_E might be forced to be equal). Brian O'Meara would consider allowing this, while cautioning you about the risks of too many parameters.
}

\value{
\code{generateMiSSEGreedyCombinations} returns a data.frame to pass to MiSSEGreedy().
}

\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models
of trait-dependent speciation and extinction. Syst. Biol. 65:583-601.

FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating trait-dependent speciation
and extinction rates from incompletely resolved phylogenies. Syst. Biol. 58:595-611.

Herrera-Alsina, L., P. van Els, and R.S. Etienne. 2018. Detecting the dependence of diversification on multiples traits from phylogenetic trees and trait data. Systematic Biology, 68:317-328.

Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary characters effect on
speciation and extinction. Syst. Biol. 56:701-710.

Magallon S. and Sanderson M.J. 2001. Absolute diversification rates in angiosperm clades. Evolution 55:1762-1780.

Nee S., May R.M., and Harvey P.H. 1994. The reconstructed evolutionary process. Philos.
Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
}

\author{Brian C. O'Meara}
\keyword{models}
