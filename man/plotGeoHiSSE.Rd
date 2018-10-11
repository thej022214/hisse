\name{plot.geohisse.states}
\alias{plot.geohisse.states}
\title{Plotting function geohisse.states objects}
\description{A plotting function for visualizing changes in states and rates over a phylogeny}
\usage{
\method{plot}{geohisse.states}(x, rate.param="net.div", do.observed.only=TRUE, rate.colors=NULL, state.colors=NULL, edge.width=5, width.factor=0.5, outline=FALSE, outline.color="grey", mar = c(0.1,0.1,0.1,0.1), type="fan", rate.range=NULL, show.tip.label=FALSE, swap.underscore=TRUE, fsize=1.0, lims.percentage.correction=0.001, legend=TRUE, legend.position=c(0, 0.2, 0, 0.2), legend.cex=0.4, legend.kernel="auto", legend.bg="cornsilk3")
}       
\arguments{
\item{x}{a \code{hisse.states} object or a list of such objects.}
\item{rate.param}{indicates the type of rates to plot. Options include: \dQuote{turnover}, \dQuote{net.div}, \dQuote{speciation},  \dQuote{extinction},  \dQuote{extinction.fraction}.}
\item{do.observed.only}{a logical indicating whether just the states should be plotted; for now, only TRUE works.}
\item{rate.colors}{user specified colors to be used for coloring rates.}
\item{state.colors}{user specified colors to be used for coloring
  states. A vector with the color for state 1 (endemic area '0'), state 2 (endemic area '1') and state 0 (widespread area '01'). Always in this order!}
\item{edge.width}{the width of the branches of the phylogeny.} 
\item{width.factor}{the factor multiplying the edge.width to get the
  width of the branch for the states tree. Needs to be a numeric value
  between 0 and 1.}
\item{type}{a character string specifying the type of phylogeny to be drawn. Options are "fan" (default) or "phylogram".}
\item{rate.range}{an optional two element vector. If present, specifies the range of rate values to use for plotting.}
\item{show.tip.label}{a logical indicating whether tip names should be included.}
\item{fsize}{sets the font size for the tip labels.}
\item{lims.percentage.correction}{deals with cases where the limits are slightly smaller than the values due to imprecision issues.}
\item{legend}{indicates if the legend is to be plotted. TRUE or FALSE.}
\item{legend.position}{the coordinates for placing the legend.}
\item{legend.cex}{the text size inside the legend.}
\item{legend.kernel}{lets you specify the way the density plot or histogram is made for rates. \dQuote{auto}
chooses what we think is the best option given your data, \dQuote{hist} makes a histogram, \dQuote{rectangular}, \dQuote{gaussian}, and others make a density plot. See \code{?density} for all non-\dQuote{hist} options.}
\item{legend.bg}{sets the color for the legend background.}
\item{mar}{sets the margins for the plot. See more information in 'help(par)'.}
\item{outline}{whether to plot an outline on the branches of the tree. Choose the color of the outline using the parameter outline.color.}
\item{outline.color}{the color for the outline. Defaults to "black".}
\item{swap.underscore}{whether to substitute all "_" with " " on the labels for the tips.}
}
\details{
This function is very similar to the \code{hisse::plot.hisse.states}
function. See more details in help page for \code{hisse::plot.hisse.states} function.
}
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 
}
\author{Brian O'Meara and Daniel Caetano}
\keyword{plotting}
