\name{Compute Model Weights}
\alias{GetAICWeights}
\title{Compute model weights}
\description{Computes the Akaike model weights for a list of HiSSE, MiSSE,
  and/or GeoHiSSE model fits.}
\usage{
GetAICWeights(hisse.results, criterion="AIC")
}
\arguments{
\item{hisse.results}{A list of models (such as from MiSSEGreedy or just putting individual hisse runs in a list)}
\item{criterion}{Which criterion to use, AIC or AICc (with correction for sample size)}
}
\details{
Function computes the model weight from their AIC values using the
formula:


\dQuote{delta <- mod.AIC - min( mod.AIC )}

\dQuote{AICw <- exp( -0.5 * delta) / sum( exp( -0.5 * delta) )}


Function will return vector of weights

}
\author{Daniel Caetano}
\keyword{AICw}
