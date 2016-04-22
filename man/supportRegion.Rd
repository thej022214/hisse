\name{SupportRegion}
\alias{SupportRegion}
\title{Adaptive Sampling of the Likelihood Surface}
\description{Adaptively samples points for each parameter to obtain an estimate of the confidence intervals.}
\usage{
SupportRegion(hisse.obj, n.points=1000, scale.int=0.1, desired.delta=2, 
output.type="turnover", hidden.states=TRUE, condition.on.survival=TRUE, 
root.type="madfitz", root.p=NULL, verbose=TRUE)
}       
\arguments{
\item{hisse.obj}{an object of class \code{hisse.fit} that contains the MLE from a model run.}
\item{n.points}{indicates the number of points to sample.}
\item{scale.int}{the scaling multiplier that defines interval to randomly sample. By default the value is set to 0.1, meaning that values are drawn at random along an interval that encompasses 10 percent above and below the MLE.}
\item{desired.delta}{defines the number lnL units away from the MLE to include. By default the value is set to 2.} 
\item{output.type}{indicates whether the rates should be printed onscreen as the optimized variables, \dQuote{turnover}, transformed to reflect net diversification, \dQuote{net.div}, or transformed to reflect \eqn{\lambda}{lambda} and \eqn{\mu}{mu}, \dQuote{raw}.}
\item{hidden.states}{a logical indicating whether the model includes a hidden state. The default is \code{TRUE}.} 
\item{condition.on.survival}{a logical indicating whether the likelihood was conditioned on the survival of two lineages and the speciation event subtending them (Nee et al. 1994). The default is \code{TRUE}.}
\item{root.type}{indicates whether root prior assumption should based the procedure described by FitzJohn et al. 2009, \dQuote{madfitz}, assumed equal, \dQuote{equal}, or set to user, \dQuote{user}.}
\item{root.p}{a vector indicating fixed root state probabilities. The default is \code{NULL}.}
\item{verbose}{a logical indicating whether progress should be printed to the screen. The default is \code{TRUE}.}
}
\details{
This function provides a means for sampling the likelihood surface quickly to estimate confidence intervals that reflect the uncertainty in the MLE. The function starts with the MLE from the hisse run. It then uses a scaling multiplier to generate an interval by which to randomly alter each parameter. However, the algorithm was designed to \dQuote{feel} the boundaries of the random search. In other words, when the algorithm begins to sample the hinterlands of the surface, it will know to restrict the boundary to allow sampling of more reasonable values based on the currently sampled set. The goal of this sampling process is to find points within some desired distance from the MLE; by default we assume this distance is 2 lnLik. The confidence interval can be estimated directly from these points. The full set of points tried are also provided and can be used to generate contour plots (though, it is not entirely straightforward to do so -- but certainly doable).
}
\value{
\code{SupportRegion} returns an object of class \code{hisse.support}. This is a list with elements:
\item{$ci}{the sampled confidence interval.}
\item{$points.within.region}{the sampled points that within 2lnL units from the MLE.}
\item{$all.points}{all points sampled by the adaptive sampler.}
}
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. In press. 

FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Syst. Biol. 58:595-611.

Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary characters effect on speciation and extinction. Syst. Biol. 56:701-710.

Nee S., May R.M., and Harvey P.H. 1994. The reconstructed evolutionary process. Philos. Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
}
\author{Jeremy M. Beaulieu}
\keyword{models}
