\name{BisseToHisse}
\alias{BisseToHisse}
\title{BisseToHisse}
\description{A simple utility function for organizing parameters from diversitree to be in the right format for \code{marginRecon} function.}
\usage{
BisseToHisse(lambda, mu, q01, q10)
}       
\arguments{
\item{lambda}{a vector containing the parameter estimates for lambda0 and lambda1 from \code{diversitree}.}
\item{mu}{a vector containing the parameter estimates for mu0 and mu1 from \code{diversitree}.}
\item{q01}{the estimate of transition rate from state 0 to state 1 from \code{diversitree}.}
\item{q10}{the estimate of transition rate from state 1 to state 0 from \code{diversitree}.} 
}
\value{
A simple vector providing the proper format for BiSSE parameters and defaults in HiSSE for use in the \code{marginRecon} function.
}
\examples{
#BiSSE ancestral reconstruction
#bisse.pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
#hisse.par.vec <- BisseToHisse(lambda=c(0.1,0.2), mu=c(0.03,0.03), q01=0.01, q10=0.01)
}
\author{Jeremy M. Beaulieu}
\keyword{models}
