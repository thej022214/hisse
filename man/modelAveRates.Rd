\name{Model averaged rates}
\alias{GetModelAveTipRates}
\alias{GetModelAveNodeRates}
\title{Model average rates at tips and nodes}
\description{Summarizes reconstructions from a single model or a set of models and returns model averaged rates.}
\usage{
GetModelAveTipRates(x)
GetModelAveNodeRates(x)
}       
\arguments{
\item{x}{a \code{hisse.states} object or a list of such objects.}
}
\details{
Provides a data frame model-averaged rates for all possible configurations of the model parameters (i.e., turnover, net.div, speciation, extinction, or extinction.fraction), either for all tips or for all nodes.
As with the plotting function, if you give a single hisse.state object, it uses that, and the rates account for uncertainty in the marginal probabilities of the states; if you give it a list of them, it will model-average the results (it assumes the trees are the same) in addition to accounting for uncertainty in the marginal probabilities at tips and nodes.
}
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 
}
\author{Jeremy M. Beaulieu}
\keyword{plotting}
