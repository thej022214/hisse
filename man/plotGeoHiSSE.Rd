\name{plot.geohisse.states}
\alias{plot.geohisse.states}
\title{Plotting function geohisse.states objects}
\description{A plotting function for visualizing changes in states and rates over a phylogeny}
\usage{
\method{plot}{geohisse.states}(x, rate.param = "net.div", type = "fan",
show.tip.label = TRUE, fsize = 1.0, legend = TRUE, ...)
}       
\arguments{
\item{x}{a \code{hisse.states} object or a list of such objects.}
\item{rate.param}{indicates the type of rates to plot. Options include: \dQuote{turnover}, \dQuote{net.div}, \dQuote{speciation},  \dQuote{extinction},  \dQuote{extinction.fraction}.}
\item{type}{a character string specifying the type of phylogeny to be drawn. Options are "fan" (default) or "phylogram".}
\item{show.tip.label}{a logical indicating whether tip names should be included. Default is TRUE}
\item{fsize}{sets the font size for the tip labels.}
\item{legend}{indicates if the legend is to be plotted. TRUE or FALSE.}
\item{...}{Additional parameters to control the plot. See \dQuote{Details}.}
}
\details{
  Additional parameters can be defined using \dQuote{...}:
  
  \dQuote{do.observed.only} is a logical indicating whether just the
  states should be plotted; for now, only TRUE
  works.
  \dQuote{rate.colors} are user specified colors to be used for
  coloring rates.
  \dQuote{state.colors} are user specified colors to be
  used for coloring states. A vector with the color for state 1 (endemic
  area '0'), state 2 (endemic area '1') and state 0 (widespread area
  '01'). Always in this order!
  \dQuote{edge.width} is the width of the
  branches of the phylogeny.
  \dQuote{width.factor} is the factor
  multiplying the \dQuote{edge.width} to get the width of the branch
  for the states tree. Needs to be a numeric value between 0 and 1.
  \dQuote{rate.range} is an optional two element
  vector. If present, specifies the range of rate values to use for
  plotting.
  \dQuote{lims.percentage.correction} deals with cases where
  the limits are slightly smaller than the values due to imprecision
  issues.
  \dQuote{legend.position} are the coordinates for placing the
  legend.
  \dQuote{legend.cex} is the text size inside the
  legend.
  \dQuote{legend.kernel} lets you specify the way the density
  plot or histogram is made for rates. A value of \dQuote{auto} chooses
  what we think is the best option given your data, \dQuote{hist} makes
  a histogram, \dQuote{rectangular}, \dQuote{gaussian}, and others make
  a density plot. See \code{?density} for all non-\dQuote{hist}
  options.
  \dQuote{legend.bg} sets the color for the legend
  background.
  \dQuote{mar} sets the margins for the plot. See more
  information in 'help(par)'.
  \dQuote{outline} is whether to plot an
  outline on the branches of the tree. Choose the color of the outline
  using the parameter outline.color.
  \dQuote{outline.color} is the color
  for the outline. Defaults to "black".
  \dQuote{swap.underscore} sets
  whether to substitute all "_" with " " on the labels for the tips.
  
  This function is very similar to the \code{hisse::plot.hisse.states}
function. See more details in help page for \code{hisse::plot.hisse.states} function.
}
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 
}
\author{Brian O'Meara and Daniel Caetano}
\keyword{plotting}
