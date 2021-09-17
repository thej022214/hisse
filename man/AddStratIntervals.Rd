\name{AddStratIntervals}
\alias{AddStratIntervals}
\title{Adds sampled stratigraphic intervals to a complete phylogeny}
\description{Takes a sampled set of stratigraphic intervals from \code{GeteStratigraphicIntervals} and adds them as segments to a plot of the phylogeny from which they were sampled.}
\usage{
AddStratIntervals(phy, f, ...)
}       
\arguments{
\item{phy}{a complete phylogeny that includes both extant and extinct lineages.} 
\item{f}{the table of sampled lineages under the fossilized birth-death process obtained from \code{GetFossils}.}
\item{...}{Additional parameters to control the points. See \code{segments}.}
}
\references{
Vascancelos, T, B.C. O'Meara, and J.M. Beaulieu. In prep.
}
\author{Jeremy M. Beaulieu}
\keyword{plotting}

