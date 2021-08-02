\name{MiSSENet}
\alias{MiSSENet}
\title{A function for catching bad model fits from MiSSEGreedy}
\description{A function for assessing whether models optimized well, identifying bad fits, and rerunning them if needed.
}
\usage{
MiSSENet(misse.list, n.tries=2, remove.bad=TRUE, dont.rerun=FALSE, save.file=NULL, 
n.cores=1, sann=TRUE, sann.its=5000, sann.temp=5230, bounded.search=TRUE, 
starting.vals=NULL, turnover.upper=10000, eps.upper=3, trans.upper=100, restart.obj=NULL)
}
\arguments{
\item{misse.list}{a \code{misse} list of model fits obtain from \code{MiSSEGreedy}.}
\item{n.tries}{maximum number of retries for a given model.}
\item{remove.bad}{a logical indicating whether models identified as poorly optimized (even after attempting to be refit) should be removed from \code{misse.list}.}
\item{dont.rerun}{a logical indicating whether models identified as poorly optimized should be run. The default is \code{FALSE}, meaning they will be rerun. Note that if \code{dont.rerun=TRUE} and \code{remove.bad=TRUE}, the bad models will be removed from the input list; when \code{remove.bad=FALSE} only a vector of model numbers are provided.}
\item{save.file}{file to use to save the full model fits before removing the poorly optimized ones.}
\item{n.cores}{how many cores to run this on in parallel.}
\item{sann}{a logical indicating whether a two-step optimization
  procedure is to be used. The first includes a simulate annealing
  approach, with the second involving a refinement using
  \code{subplex}. The default is \code{FALSE}.}
\item{sann.its}{a numeric indicating the number of times the simulated
  annealing algorithm should call the objective function.}
\item{sann.temp}{the starting temperature for the simulated annealing. Higher temperatures results in the chain sampling a much wider space initially. The default of 5320 is based on the default of the GenSA package. For larger trees setting this value higher in conjunction with more sann.its can drastically improve performance.}
\item{bounded.search}{a logical indicating whether or not bounds should
  be enforced during optimization. The default is \code{TRUE}.}
\item{starting.vals}{a numeric vector of length 3 with starting values for the
  model. Position [1] sets turnover, [2] sets extinction fraction, and
  [3] transition rates between distinct diversification rates.}
\item{turnover.upper}{sets the upper bound for the turnover parameters.}
\item{eps.upper}{sets the upper bound for the eps parameters.}
\item{trans.upper}{sets the upper bound for the transition rate parameters.}
\item{restart.obj}{an object of class  that contains everything to restart an optimization.}
}
\details{
This function is used to triage poorly optimized models after a \code{MiSSEGreedy} run. 
It is normally invoked within \code{MiSSEGreedy}, but it can also be used as a standalone function, 
to simply identify poorly identify models and/or rerun them.
}
\references{
Vascancelos, T, B.C. O'Meara, and J.M. Beaulieu. In prep.
}
\author{Jeremy M. Beaulieu}
\keyword{models}
