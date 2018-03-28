\name{Model averaged equilibrium frequencies}
\alias{GetModelAveEqFreqs}
\title{Model average equlibrium frequencies}
\description{Calculates the expected equilibrium frequencies under various geographic models from a single model or a set of models and returns model averaged frequencies.}
\usage{
GetModelAveEqFreqs(x, max.time)
}       
\arguments{
\item{x}{a \code{higeosse.fit} object, or a list of such objects. 
A list of model can only include models of one type.}
\item{max.time}{a numeric value that indicates the duration of time
when calculating equilibrium frequencies.}
}
\details{
Provides a data frame model-averaged equilibrium frequencies across
a set of models. If you give a single higeosse.fir object,
it uses that; if you give it a list of them, it will
model-average the results.
}
\references{
Caetano, D.S., B.C. O'Meara, and J.M. Beaulieu. In prep.
}
\author{Jeremy M. Beaulieu}
\keyword{plotting}
