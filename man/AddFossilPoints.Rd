\name{AddFossilPoints}
\alias{AddFossilPoints}
\title{Adds sampled fossils to a complete phylogeny}
\description{Takes a sampled set of m and k fossils from \code{GetFossils} and them as points to a plot of the phylogeny from which they were sampled.}
\usage{
AddFossilPoints(phy, f, ...)
}       
\arguments{
\item{phy}{a complete phylogeny that includes both extant and extinct lineages.} 
\item{f}{the table of sampled lineages under the fossilized birth-death process obtained from \code{GetFossils}.}
\item{...}{Additional parameters to control the points. See \code{points}.}
}
\references{
Vascancelos, T, B.C. O'Meara, and J.M. Beaulieu. In prep.
}
\author{Jeremy M. Beaulieu}
\keyword{plotting}


