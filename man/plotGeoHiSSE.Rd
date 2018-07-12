\name{plot.geohisse.states}
\alias{plot.geohisse.states}
\title{Plotting function geohisse.states objects}
\description{A plotting function for visualizing changes in states and rates over a phylogeny}
\usage{
\method{plot}{geohisse.states}(x, rate.param, do.observed.only=TRUE, rate.colors=NULL, 
state.colors=NULL, edge.width.rate=5, edge.width.state=2, type="fan", 
rate.range=NULL, show.tip.label=TRUE, fsize=1.0, lims.percentage.correction=0.001,
legend="tips", legend.position=c(0, 0.2, 0, 0.2), legend.cex=0.4, 
legend.kernel.rates="auto", legend.kernel.states="auto", 
legend.bg="cornsilk3",...)
}       
\arguments{
\item{x}{a \code{hisse.states} object or a list of such objects.}
\item{rate.param}{indicates the type of rates to plot. Options include: \dQuote{turnover}, \dQuote{net.div}, \dQuote{speciation},  \dQuote{extinction},  \dQuote{extinction.fraction}.}
\item{do.observed.only}{a logical indicating whether just the states should be plotted; for now, only TRUE works.}
\item{rate.colors}{user specified colors to be used for coloring rates.}
\item{state.colors}{user specified colors to be used for coloring states.} 
\item{edge.width.rate}{the width of the rate lines.} 
\item{edge.width.state}{the width of the state lines, which should be set so that it is narrower than the edge width for rates.}
\item{type}{a character string specifying the type of phylogeny to be drawn. Options are "fan" (default) and "phylogram". See
  \code{ape::plot.phylo} function for more options.}
\item{rate.range}{an optional two element vector. If present, specifies the range of rate values to use for plotting.}
\item{show.tip.label}{a logical indicating whether tip names should be included.}
\item{fsize}{sets the font size for the tip labels.}
\item{lims.percentage.correction}{deals with cases where the limits are slightly smaller than the values due to imprecision issues.}
\item{legend}{indicates the type of legend. Options include: \dQuote{none}, \dQuote{traditional}, \dQuote{tips}, \dQuote{internal}, \dQuote{all}.}
\item{legend.position}{the coordinates for placing the legend.}
\item{legend.cex}{the text size inside the legend.}
\item{legend.kernel.rates}{for legend=tips, internal, or all, lets you specify the way the density plot or histogram is made for rates. \dQuote{auto}
chooses what we think is the best option given your data, \dQuote{hist} makes a histogram, \dQuote{rectangular}, \dQuote{gaussian}, and others make a density plot. See \code{?density} for all non-\dQuote{hist} options.}
\item{legend.kernel.states}{as above, for states.}
\item{legend.bg}{sets the color for the legend background.}
\item{\dots}{further arguments to be passed to \dQuote{plot} or \dQuote{plot.hisse.states}.}
}
\details{
This function is very similar to the \code{hisse::plot.hisse.states}
function. See more details in help page for \code{hisse::plot.hisse.states} function.
}
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 
}
\author{Brian O'Meara}
\keyword{plotting}
