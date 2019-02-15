\name{Model averaged equilibrium frequencies}
\alias{GetModelAveEqFreqs}
\title{Model average equlibrium frequencies}
\description{Calculates the expected equilibrium frequencies under various hidden SSE models from a single model or a set of models and returns model averaged frequencies.}
\usage{
GetModelAveEqFreqs(x, max.time, model.type="hisse", get.rates=FALSE, 
rate.type="turnover", get.all.states=FALSE)
}       
\arguments{
\item{x}{a \code{
.fit} object, or a list of such objects. 
A list of model can only include models of one type.}
\item{max.time}{a numeric value that indicates the duration of time
when calculating equilibrium frequencies. 
(Not the maximum time that the function is allowed to run!)}
\item{model.type}{indicates which set of models to calculate equilibrium frequencies for. Options are "hisse", "geohisse", or "muhisse"}
\item{get.rates}{indicates whether you want model-averaged rates based on equilibrium. Default is FALSE, which means just the equilibrium frequencies are returned.}
\item{rate.type}{if get.rates=TRUE, this specifies which type of rates to be returned. Default is "turnover", but options include "net.div" and "speciation".}
\item{get.all.states}{indicates that you would like the equilibrium frequencies to be output for all observed and hidden state combinations. Note that under this option, the frequencies cannot be averaged across models.} 
}
\details{
Provides a data frame model-averaged equilibrium frequencies across a set of models. If, for example, you give it a single geohisse.fit object, it uses that; if you give it a list of them, it will model-average the results. The purpose of this approach is dealing with mixed signals among hidden states with respect to the observed states. We discuss this in detail in Caetano et al (2018).

Note, that this function does not work for model fits obtained from fGeoHiSSE. This will provided in the next version.  
}
\references{
Caetano, D.S., B.C. O'Meara, and J.M. Beaulieu. 2018. Hidden state models improve state-dependent diversification approaches, including biogeographic models. Evolution, 72:2308-2324.
}
\author{Jeremy M. Beaulieu}
\keyword{plotting}
