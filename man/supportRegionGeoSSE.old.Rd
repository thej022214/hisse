\name{SupportRegionGeoSSE.old}
\alias{SupportRegionGeoSSE.old}
\title{Adaptive Sampling of the Likelihood Surface under the original GeoHiSSE model}
\description{Adaptively samples points for each parameter to obtain an estimate of the confidence intervals.}
\usage{
SupportRegionGeoSSE.old(geohisse.old.obj, n.points=1000, scale.int=0.1, desired.delta=2, 
min.number.points=10, verbose=TRUE)
}       
\arguments{
\item{geohisse.old.obj}{an object of class \code{geohisse.old.fit} that contains the MLE from a model run.}
\item{n.points}{indicates the number of points to sample.}
\item{scale.int}{the scaling multiplier that defines interval to randomly sample. By default the value is set to 0.1, meaning that values are drawn at random along an interval that encompasses 10 percent above and below the MLE.}
\item{desired.delta}{defines the number lnL units away from the MLE to include. By default the value is set to 2.} 
\item{min.number.points}{sets the minimum number of points that can be returned. By default the value is set to 10.}
\item{verbose}{a logical indicating whether progress should be printed to the screen. The default is \code{TRUE}.}
}
\details{
This function provides a means for sampling the likelihood surface quickly to estimate confidence intervals that reflect the uncertainty in the MLE. The function starts with the MLE from the hisse run. It then uses a scaling multiplier to generate an interval by which to randomly alter each parameter. However, the algorithm was designed to \dQuote{feel} the boundaries of the random search. In other words, when the algorithm begins to sample the hinterlands of the surface, it will know to restrict the boundary to allow sampling of more reasonable values based on the currently sampled set. The goal of this sampling process is to find points within some desired distance from the MLE; by default we assume this distance is 2 lnLik. The confidence interval can be estimated directly from these points. The full set of points tried are also provided and can be used to generate contour plots (though, it is not entirely straightforward to do so -- but certainly doable).

Note the scale.int option. This roughly sets the variance for sampling points. If this seems to take a long while to find enough points within the desired likelihood region consider reducing scale.int to either 0.05 or, in some cases, 0.01.
}
\value{
\code{SupportRegionGeoSSE.old} returns an object of class \code{geohisse.old.support}. This is a list with elements:
\item{$ci}{the sampled confidence interval.}
\item{$points.within.region}{the sampled points that within 2lnL units from the MLE.}
\item{$all.points}{all points sampled by the adaptive sampler.}
}
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 

FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Syst. Biol. 58:595-611.

Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary characters effect on speciation and extinction. Syst. Biol. 56:701-710.

Nee S., May R.M., and Harvey P.H. 1994. The reconstructed evolutionary process. Philos. Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
}
\author{Jeremy M. Beaulieu}
\keyword{models}
