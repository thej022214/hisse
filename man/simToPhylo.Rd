\name{SimToPhylo}
\alias{SimToPhylo}
\title{Convert simulated result to a tree}
\description{Converts the $results element of the return to a phylo object}
\usage{
SimToPhylo(results, include.extinct=FALSE, drop.stem=TRUE)
}       
\arguments{
\item{results}{dataframe of results from SimulateHisse.}
\item{include.extinct}{include (TRUE) or delete (FALSE) extinct terminal taxa.}
\item{drop.stem}{include the lineage leading to the first saved speciation even (FALSE) or cull it (TRUE).}
}
\details{
This takes the \code{$results} object from a \code{SimulateHisse} run and converts it to an ape \code{phylo} object. If there are no taxa on the final tree, it returns NA; if there is one taxon, it returns a one taxon tree. This is behavior different from diversitree's tree simulators, which returns NULL in both the zero and one taxon case. Extinct taxa can be pruned or not. For simulations starting with one taxon, and/or for simulations with some extinction, the final tree can have a time when there is a single lineage before it radiates into the crown group. This stem can be included or not. The tip states are stored in \code{$tip.state} in the returned tree.
}
\value{
\item{phylo}{a phylo object in cladewise order.}
}
\examples{
\dontrun{
simulated.result <- SimulateHisse(c(.3, .1), c(.1, 0), 
matrix(c(NA, 0.2, .3, NA), nrow=2), max.taxa=35, x0=1)	
par(mfcol=c(1,2))
plot(SimToPhylo(simulated.result$results, include.extinct=TRUE))
plot(SimToPhylo(simulated.result$results, include.extinct=FALSE))
} 	
}
\references{
Beaulieu, J.M, and B.C. O'Meara. In revision. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. In revision. 
}
\author{Brian O'Meara}
\keyword{utility}
