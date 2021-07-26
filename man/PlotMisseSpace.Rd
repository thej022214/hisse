\name{PlotMisseSpace}
\alias{PlotMisseSpace}
\title{Plotting function for MiSSEGreedy}
\description{A plotting function for visualizing the model space explored in MiSSEGreedy}
\usage{
\method{PlotMisseSpace}(x, possible.combos, ...)
}       
\arguments{
\item{x}{a \code{misse.states} list of such objects.}
\item{possible.combos}{data.frame of parameter combinations to try. See 'Details'.}
\item{...}{Additional parameters to control the igraph plot. See \dQuote{Details}.}
}
\details{
  Additional parameters can be defined using \dQuote{...}:
  
  \dQuote{do.observed.only} is a logical indicating whether just the
  states should be plotted; for now, only TRUE
  works.
  \dQuote{rate.colors} are user specified colors to be used for
  coloring rates.
  \dQuote{state.colors} are user specified colors to be
  used for coloring states. A vector with the color for states 00, 01, 10, and 11 (in this order).
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
Vascancelos, T, B.C. O'Meara, and J.M. Beaulieu. In prep.
\author{Jeremy M. Beaulieu}
\keyword{plotting}
