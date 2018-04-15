\name{TranslateParsMakerGeoHiSSE}
\alias{TranslateParsMakerGeoHiSSE}
\title{Key to translate GeoHiSSE into ClaSSE}
\description{Simulate key of parnames to translate between GeoHiSSE and ClaSSE models}
\usage{
TranslateParsMakerGeoHiSSE(k, add.extinction, add.jumps)}       
\arguments{
\item{k}{number of hidden states in the model. Minimum is 0 for a model
  without hidden states (i.e., GeoSSE).}
\item{add.extinction}{if the extinction parameter should be separated from extirpation.}
\item{add.jumps}{if jumps between endemic areas should be added.}
}
\details{
Function will only return the parameters that are relevant for a
GeoHiSSE model. The ClaSSE model will have many more parameters. The
extra parameters of the ClaSSE model need to be all set to 0. The number
of parameters of the correspondent ClaSSE model given the parameter
\code{k} is given by:

(k+1)*3

}
\value{
  Function returns a matrix of characters names. Column 'classe.pars'
  have the parameter names for a ClaSSE model and column 'geohisse.pars'
  have the parameter names for the matching GeoHiSSE model.
      }
\author{Daniel Caetano}
\keyword{simulation, GeoSSE, GeoHiSSE}
