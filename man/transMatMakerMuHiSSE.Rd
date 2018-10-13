\name{TransMatMakerMuHiSSE}
\alias{TransMatMakerMuHiSSE}
\title{Transition Rate matrix generator for the MuHiSSE model}
\description{Generates and manipulates the index of the transition rate parameters to be optimized.}
\usage{
TransMatMakerMuHiSSE(hidden.traits=0, make.null=FALSE, include.diagonals=FALSE,
make.special.for.teo=FALSE)
}       
\arguments{
\item{hidden.traits}{a numeric value with the number of hidden states in
  the model. The cannonical MuSSE model has no hidden states, so
  \code{'hidden.states=0'}. The default value is 0.} 
\item{make.null}{Sets the transition matrix to the null model such that (1)
  transition rates are the same across all hidden states and (2) the transition rates
  between hidden states (such as 00A <-> 01A) are the same.}
\item{include.diagonals}{allows for dual transitions, e.g., 00->11 and 11->00 transitions.}
\item{make.special.for.teo}{Testing proof of concept for Teo.}
}
\details{
  Outputs the full index of the rate parameters that are to be optimized.

  The intention is that a user might want to see how the matrix is designed prior to an analysis and perhaps drop a few parameters beforehand due to some hypothesis that she might have. The resulting matrix is to be plugged directly into \code{MuHiSSE}.
}
\value{
Returns a rate matrix index.
}
