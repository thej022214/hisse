\name{plot.hisse.states}
\alias{plot.hisse.states}
\title{Plotting function hisse.states objects}
\description{A plotting function for visualizing changes in states and rates over a phylogeny}
\usage{
\method{plot}{hisse.states}(x, rate.param, do.observed.only=TRUE, rate.colors=NULL, 
state.colors=NULL, edge.width=5, width.factor=0.5, type="fan", 
mar=c(0.1,0.1,0.1,0.1), outline=FALSE, outline.color="black",
rate.range=NULL, show.tip.label=TRUE, swap.underscore=TRUE, 
fsize=1.0, lims.percentage.correction=0.001,
legend="tips", legend.position=c(0, 0.2, 0, 0.2), legend.cex=0.4, 
legend.kernel.rates="auto", legend.kernel.states="auto", 
legend.bg="cornsilk3")
}       
\arguments{
\item{x}{a \code{hisse.states} object or a list of such objects.}
\item{rate.param}{indicates the type of rates to plot. Options include: \dQuote{turnover}, \dQuote{net.div}, \dQuote{speciation},  \dQuote{extinction},  \dQuote{extinction.fraction}.}
\item{do.observed.only}{a logical indicating whether just the states should be plotted; for now, only TRUE works.}
\item{rate.colors}{user specified colors to be used for coloring rates.}
\item{state.colors}{user specified colors to be used for coloring states.} 
\item{edge.width}{the width of the branches of the phylogeny.} 
\item{width.factor}{the factor multiplying the edge.width to get the
  width of the branch for the states tree. Needs to be a numeric value
  between 0 and 1.}
\item{type}{a character string specifying the type of phylogeny to be drawn. Options are "fan" (default) or "phylogram".}
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
\item{mar}{sets the margins for the plot. See more information in 'help(par)'.}
\item{outline}{whether to plot an outline on the branches of the tree. Choose the color of the outline using the parameter outline.color.}
\item{outline.color}{the color for the outline. Defaults to "black".}
\item{swap.underscore}{whether to substitute all "_" with " " on the labels for the tips.}
}
\details{
Provides phylogeny that shows a heat map of the diversification rate parameter you specify (which could be turnover, net.div, speciation, extinction, or extinction.fraction).
The discrete state reconstruction appears as lines on top of the heat map. If you give a single hisse.state object, it uses that; if you give it a list of them, it will 
model-average the results (it assumes the trees are the same). Colors can be specified by sending a vector of colors to rate.colors or state.colors (the defaults
are red to blue for rate and white to black for state). You can specify two or more colors: c("red", "gray", "blue") for example. By default the visualization uses the
minimum rate on the tree for the minimum color, and the maximum rate for the maximum color, but you may want to use the same color scale across models, even
if some of them have a smaller range than others. To do this, pass a vector with the minimum and maximum rate across all models to the visualization for all models
and they will use the same scale. There are many options for adding a legend. A traditional legend showing what values a color corresponds to is \dQuote{traditional}, like
what plotSimmap will show in phytools. However, we can also use the legend to show the distribution of values, rather than just a key to color. \dQuote{tips} shows a density
plot of states or rates at tips, \dQuote{internal} a distribution at internal nodes, and \dQuote{all} at all nodes in the tree. For the density or histogram plots, you can
let the package pick the best visualization or choose yourself whether to use a histogram or density plot, and if the latter, what kernel you want. The legend can be
moved around the overall tree plot by using \dQuote{legend.position}: this is a vector that specifies the \dQuote{fig} argument to \dQuote{par}: c(x1, x2, y1, y2), where
the values are the starting and ending positions as a fraction of the overall plot. By default, the legend starts at the lower left corner and continues up 20% into 
the rest of the plot (c(0, 0.2, 0, 0.2)): by changing values, you can make the legend larger or smaller and change its position. The heatmap code is modified slightly
from Liam Revell's phytools package.
}
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 
}
\author{Brian O'Meara and Daniel Caetano}
\keyword{plotting}
