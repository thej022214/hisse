\name{PlotMisseSpace}
\alias{PlotMisseSpace}
\title{Plotting function for MiSSEGreedy}
\description{A plotting function for visualizing the model space explored in MiSSEGreedy. 
Arguments in \code{...} are passed to \code{plot.igraph}.
}
\usage{
PlotMisseSpace(x, possible.combos=NULL, arrows.by.weight=FALSE, ...)
}       
\arguments{
\item{x}{a \code{misse} list of model fits obtain from \code{MiSSEGreedy}.}
\item{possible.combos}{data.frame of parameter combinations to try. See 'Details'.}
\item{arrows.by.weight}{a logical indicating whether arrow direction between adjacent model reflects which model has higher support. The default is \code{FALSE} which means arrows show direction of model complexity.}
\item{...}{Additional parameters to control the igraph plot. See \code{plot.igraph}.}
}
\details{
If the input x is a list of \code{misse} fits from MiSSEGreedy then size of the vertices are in proportion 
of their AIC weights. If only the \dQuote{possible.combos} is input, and \code{x=NULL}, then the model space 
will be plotted but with each vertex having the same size.
}
\references{
Vascancelos, T, B.C. O'Meara, and J.M. Beaulieu. In prep.
}
\author{Jeremy M. Beaulieu}
\keyword{plotting}
