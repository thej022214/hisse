\name{PlotMisseSpace}
\alias{PlotMisseSpace}
\title{Plotting function for MiSSEGreedy}
\description{A plotting function for visualizing the model space explored in MiSSEGreedy}
\usage{
PlotMisseSpace(x, possible.combos, ...)
}       
\arguments{
\item{x}{a \code{misse} list of model fits obtain from \code{MiSSEGreedy}.}
\item{possible.combos}{data.frame of parameter combinations to try. See 'Details'.}
\item{...}{Additional parameters to control the igraph plot. See \code{igraph}.}
}
\details{
If the input x is a list of \code{misse} fits from MiSSEGreedy then size of the vertices are in proportion 
of their AIC weights. If only the \dQuote{possible.combos} is input, and \code{x=NULL}, then the model space 
will be plotted but with each vertices having the same size.
}
\references{
Vascancelos, T, B.C. O'Meara, and J.M. Beaulieu. In prep.
}
\author{Jeremy M. Beaulieu}
\keyword{plotting}
