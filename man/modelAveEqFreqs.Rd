\name{Model averaged equilibrium frequencies}
\alias{GetModelAveEqFreqs}
\title{Model average equlibrium frequencies}
\description{Calculates the expected equilibrium frequencies under various geographic models from a single model or a set of models and returns model averaged frequencies.}
\usage{
GetModelAveEqFreqs(x, max.time, model.type="hisse", get.rates=FALSE, 
rate.type="turnover")
}       
\arguments{
\item{x}{a \code{
.fit} object, or a list of such objects. 
A list of model can only include models of one type.}
\item{max.time}{a numeric value that indicates the duration of time
when calculating equilibrium frequencies.}
\item{model.type}{indicates which set of models to calculate equilibrium frequencies for. Options are "hisse" or "geohisse"}
\item{get.rates}{indicates whether you want model-averaged rates based on equilibrium. Default is FALSE, which means just the equilibrium frequencies are returned.}
\item{rate.type}{if get.rates=TRUE, this specifies which type of rates to be returned. Default is "turnover", but options include "net.div" and "speciation".}
}
\details{
Provides a data frame model-averaged equilibrium frequencies across
a set of models. If you give a single geohisse.fit object,
it uses that; if you give it a list of them, it will
model-average the results.
}
\references{
Caetano, D.S., B.C. O'Meara, and J.M. Beaulieu. In prep.
}
\author{Jeremy M. Beaulieu}
\keyword{plotting}
