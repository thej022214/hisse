\name{TransMatMakerGeoHiSSE}
\alias{TransMatMakerGeoHiSSE}
\title{Transition Rate matrix generator for the GeoHiSSE model}
\description{Generates and manipulates the index of the transition rate parameters to be optimized.}
\usage{
TransMatMakerGeoHiSSE(hidden.areas=0, make.null=FALSE, include.jumps=FALSE, 
	separate.extirpation=FALSE)
}       
\arguments{
\item{hidden.areas}{a numeric value with the number of hidden areas in
  the model. The cannonical GeoSSE model has no hidden areas, so
  \code{'hidden.areas=0'}. The default value is 0.} 
\item{make.null}{Sets the transition matrix to the null model such that (1)
  transition rates from endemic to widespread areas (A -> AB and B ->
  AB) are the same across all hidden areas and (2) the transition rates
  between hidden areas (such as 0A <-> 1A) are the same for all areas.}
\item{include.jumps}{allows for 0->1 and 1->0 transitions.}
\item{separate.extirpation}{allows for 01->0 and 01->1 transitions.}
}
\details{
  Outputs the full index of the rate parameters that are to be optimized.

  The intention is that a user might want to see how the matrix is designed prior to an analysis and perhaps drop a few parameters beforehand due to some hypothesis that she might have. The resulting matrix is to be plugged directly into \code{GeoHiSSE}.
}
\value{
Returns a rate matrix index.
}
