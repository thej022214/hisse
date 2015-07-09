\name{SimulateHisse}
\alias{SimulateHisse}
\title{Simulate under a HiSSE model}
\description{Flexible simulation function allowing checkpointing}
\usage{
\method{SimulateHisse}(turnover.rates, eps.values, transition.rates, max.taxa=Inf, max.t=Inf, max.wall.time=Inf, x0, nstart=1, checkpoint.file=NULL, checkpoint.frequency=100, checkpoint.start.object=NULL, override.safeties=FALSE)
}       
\arguments{
\item{turnover.rates}{a vector of turnover rates.}
\item{eps.values}{a vector of extinction fractions.}
\item{transition.rates}{a matrix of transition rates.}
\item{max.taxa}{have the simulation stop at \code{max.taxa} surviving taxa.}
\item{max.t}{have the simulation stop at \code{max.t} height (including stem).} 
\item{max.wall.time}{have the simulation stop at this many seconds of running.} 
\item{x0}{the starting state.}
\item{nstart}{how many taxa you want to start from. Usual values are 1 (start with single lineage) or 2 (crown group).}
\item{checkpoint.file}{if you want to save progress (so you can restart from last saved point) include the name of a file.}
\item{checkpoint.frequency}{if there is a checkpoint file, frequency (in terms of number of steps: births, deaths, or transitions) at which this is saved.}
\item{checkpoint.start.object}{if you are starting from a checkpointed.file, the object loaded from that file. Typically called \code{checkpoint.result}.}
\item{override.safeties}{simulate even if there is no limit on number of taxa, depth of tree, or wall time.}
}
\details{
Note that currently, the simulator assumes turnover.rates, eps.values, and transition rates are in the same state order

Provides phylogeny that shows a heat map of the diversification rate parameter you specify (which could be turnover, net.div, speciation, extinction, or extinction.fraction).
The discrete state reconstruction appears as lines on top of the heat map. If you give a single hisse.state object, it uses that; if you give it a list of them, it will 
model-average the results (it assumes the trees are the same). Colors can be specified by sending a vector of colors to rate.colors or state.colors (the defaults
are red to blue for rate and white to black for state). You can specify two or more colors: c("red", "gray", "blue") for example. By default the visualization uses the
minimum rate on the tree for the minimum color, and the maximum rate for the maximum color, but you may want to use the same color scale across models, even
if some of them have a smaller range than others. To do this, pass a vector with the minimum and maximum rate across all models to the visualization for all models
and they will use the same scale. There are many options for adding a legend. A traditional legend showing what values a color corresponds to is dQuote{traditional}, like
what plotSimmap will show in phytools. However, we can also use the legend to show the distribution of values, rather than just a key to color. dQuote{tips} shows a density
plot of states or rates at tips, dQuote{internal} a distribution at internal nodes, and dQuote{all} at all nodes in the tree. For the density or histogram plots, you can
let the package pick the best visualization or choose yourself whether to use a histogram or density plot, and if the latter, what kernel you want. The legend can be
moved around the overall tree plot by using dQuote{legend.position}: this is a vector that specifies the dQuote{fig} argument to dQuote{par}: c(x1, x2, y1, y2), where
the values are the starting and ending positions as a fraction of the overall plot. By default, the legend starts at the lower left corner and continues up 20% into 
the rest of the plot (c(0, 0.2, 0, 0.2)): by changing values, you can make the legend larger or smaller and change its position. The heatmap code is modified slightly
from Liam Revell's phytools package.
}
\value{
}
\references{
Beaulieu, J.M, and B.C. O'Meara. In revision. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. In revision. 
}
\author{Brian O'Meara}
\keyword{simulation}
