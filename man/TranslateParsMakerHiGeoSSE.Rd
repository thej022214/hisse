\name{TranslateParsMakerHiGeoSSE}
\alias{TranslateParsMakerHiGeoSSE}
\title{Key to translate HiGeoSSE into ClaSSE}
\description{Simulate key of parnames to translate between HiGeoSSE and ClaSSE models}
\usage{
TranslateParsMakerHiGeoSSE <- function(k)}       
\arguments{
\item{k}{number of hidden states in the model. Minimum is 0 for a model
  without hidden states (i.e., GeoSSE).}
}
\details{
Function will only return the parameters that are relevant for a
HiGeoSSE model. The ClaSSE model will have many more parameters. The
extra parameters of the ClaSSE model need to be all set to 0. The number
of parameters of the correspondent ClaSSE model given the parameter
\code{k} is given by:

(k+1)*3

}
\value{
  Function returns a matrix of characters names. Column 'classe.pars'
  have the parameter names for a ClaSSE model and column 'higeosse.pars'
  have the parameter names for the matching HiGeoSSE model.
      }
\examples{}
\references{}
\author{Daniel Caetano}
\keyword{simulation, GeoSSE, HiGeoSSE}
