\name{Compute Model Weights}
\alias{GetModelWeight}
\title{Compute model weights}
\description{Computes the Akaike model weights for a list of HiSSE or
  GeoHiSSE model fits.}
\usage{
GetModelWeight(...)
}       
\arguments{
\item{...}{Any number of model fits of class 'hisse.fit' or
  'geohisse.fit' separated by commas OR a single list object composed by
the model fits.}
}
\details{
Function computes the model weight from their AIC values using the
formula:


\dQuote{delta <- mod.AIC - min( mod.AIC )}

\dQuote{aicw <- exp( -0.5 * delta) / sum( exp( -0.5 * delta) )}


Function will return model names is these are provided in the call using
'...' or a list object.    
    
}
\author{Daniel Caetano}
\keyword{AICw}