\name{SimulateHiGeoSSE}
\alias{SimulateHiGeoSSE}
\title{Simulate under a HiGeoSSE model}
\description{Simulation function for the HiGeoSSE model}
\usage{
SimulateHiGeoSSE <- function(mod.matrix, trans.matrix, max.taxa=Inf,
max.t=Inf, max.wall.time=Inf, include.extinct=FALSE, x0="AB0",
override.safeties=FALSE, return.empty.pars=FALSE, control.empty.pars=NULL)}       
\arguments{
\item{mod.matrix}{a matrix of parameter values for the model. Number of
rows equal to the parameters of the GeoSSE model and number of
columns (above 1) equal to the number of hidden states to be
simulated. Rownames have to match the parameters of the GeoSSE
model. See 'Details'.}
\item{trans.matrix}{a squared transition matrix between
states. Dimension (nrow and ncol) need to be a multiple of 3. Order of
traits are A0, B0, AB0, A1, B1, AB1 (for a single hidden trait) and so
on. See 'Details'.}
\item{max.taxa}{have the simulation stop at \code{max.taxa} surviving taxa.}
\item{max.time}{have the simulation stop at \code{max.time} height (including stem).} 
\item{max.wall.time}{have the simulation stop at this many seconds of running.} 
\item{x0}{a character value for the starting state. Default is "AB0".}
\item{override.safeties}{simulate even if there is no limit on number of
taxa, depth of tree, or wall time.}
\item{return.empty.pars}{whether to not run the simulation and return a
list with $mod.matrix and $trans.matrix in the correct format to be used
for this function.}
\item{control.empty.pars}{control the options for
'return.empty.pars'. Need to be a list with: $k - number of hidden
states, $s - speciation rate, $x - extirpation rate, $d
- dispersal rate, $q - transition rate. If a single value is provided,
then it is repeated for all parameters. If a named vector is provided
then parameters are used accordingly to matched names (see Details) and missing parameters are set to 0.}
}
\details{
The function assumes a set of patternized names for the parameters of
the model. For example, to be used as rownames of the \code{mod.matrix}
or by the \code{control.empty.pars} options. The first set of parameters
are always:

sA0, sB0, sAB0, xA0, xB0, dA0, dB0

The set of parameters for the first hidden rate are:

sA1, sB1, sAB1, xA1, xB1, dA1, dB1

and so on, in sequential number starting from 0. For the transition rates between the hidden rates, we use:

qA01, qA10, qB01, qB10, qAB01, qAB10

The \code{return.empty.pars} option will return a list of the correct
number of parameters for the model. This makes it much easier to avoid
errors when running simulations. See the examples.

This function uses \code{diversitree::tree.classe} in order to performs simulations.

We strongly advise putting some limits to stop the run. These can be any
combination of \code{max.taxa}, \code{max.t}, and \code{max.wall.time}.

}
\value{
	This returns a list, with the following elements:
	\item{$phy}{a phylogeny ot class phylo.}
	\item{$data}{a named vector with the tip states.}
	\item{$birth.counts}{a named vector with the number of
cladogenetic events for each area (and hidden area).}
	\item{$death.counts}{a named vector with the number of local
extinctions in each of the endemic areas (hidden or not).}
	\item{$extipation.counts}{a named vector with the number of
extirpation events from the widespread distribution (hidden or not).}
	\item{$dispersion.counts}{a named vector with the number of dispersion
events from the endemic areas to the widespread area (hidden or not).}
	\item{$hidden.trans.counts}{a named vector with the number of
transition events between the hidden categores for each parameter.}
	\item{$transition.counts}{matrix of counts of transitions from
	  one state to another.}
      }      
\examples{
\dontrun{
## Get the a list with the correct parameters to run the simulation:
SimulateHiGeoSSE(return.empty.pars=TRUE, control.empty.pars=list(k=1))
} 	
}
\references{}
\author{Daniel Caetano}
\keyword{simulation, GeoSSE, HiGeoSSE}
