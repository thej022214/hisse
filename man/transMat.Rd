\name{TransMatMaker}
\alias{TransMatMaker}
\alias{ParDrop}
\alias{ParEqual}
\title{Transition Rate matrix generator}
\description{Generates and manipulates the index of the rate parameters to be optimized}
\usage{
TransMatMaker(hidden.states=FALSE)
ParDrop(rate.mat.index=NULL, drop.par=NULL)
ParEqual(rate.mat.index=NULL, eq.par=NULL)
}       
\arguments{
\item{hidden.states}{a logical indicating whether the underlying model includes hidden states. The default is \code{FALSE}.} 
\item{rate.mat.index}{A user-supplied rate matrix index to be manipulated.}
\item{drop.par}{a vector of transitions to be dropped from the model.}
\item{eq.par}{a vector of transitions pairs to be set equal.}
}
\details{
Outputs the full index of the rate parameters that are to be optimized. The intention is that a user might want to see how the matrix is designed prior to an analysis and perhaps drop a few parameters beforehand due to some hypothesis that he or she might have. The resulting matrix is to be plugged directly into \code{hisse}.
}
\value{
Returns a rate matrix index
}
