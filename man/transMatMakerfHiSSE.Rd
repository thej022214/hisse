\name{TransMatMakerfHiSSE}
\alias{TransMatMakerfHiSSE}
\title{Transition Rate matrix generator for the fHiSSE model}
\description{Generates and manipulates the index of the transition rate parameters to be optimized.}
\usage{
TransMatMakerfHiSSE(hidden.traits=0, make.null=FALSE, cat.trans.vary=FALSE)
}       
\arguments{
\item{hidden.traits}{a numeric value with the number of hidden states in
  the model. The canonical BiSSE model has no hidden states, so
  \code{'hidden.traits=0'}. The default value is 0.} 
\item{make.null}{Sets the transition matrix to the null model such that (1)
  transition rates are the same across all hidden states and (2) the transition rates
  between hidden states (such as 0A <-> 1A) are the same.}
\item{cat.trans.vary}{sets transition among hidden categories to vary.}
}
\details{
  Outputs the full index of the rate parameters that are to be optimized.

  The intention is that a user might want to see how the matrix is designed prior to an analysis and perhaps drop a few parameters beforehand due to some hypothesis that she might have. The resulting matrix is to be plugged directly into \code{hisse.new}.
}
\value{
Returns a rate matrix index.
}
