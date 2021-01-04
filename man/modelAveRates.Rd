\name{Model averaged rates}
\alias{GetModelAveRates}
\title{Model average rates at tips and nodes}
\description{Summarizes reconstructions from a single model or a set of models and returns model averaged rates.}
\usage{
GetModelAveRates(x, AIC.weights=NULL, type=c("tips", "nodes", "both"),
bound.par.matrix=cbind(c(0,-1000,0,0,0),c(1000,1000,1000,3,1000)) )
}       
\arguments{
\item{x}{a \code{hisse.states} object, a \code{hisse.states.geosse}
  object, a \code{muhisse.states} or a list of such objects. A 
  list of model can only include models of one type.}
\item{AIC.weights}{a numeric vector with length equal to the number of
  models in the list 'x'. This is the AICw for each of the models to be
  averaged. If 'AIC.weights==NULL' (the default value) then function
  will compute the AICw for the models in the set from the $AIC value
  for each model.}
\item{type}{one of "tips", "nodes", or "both" (the default). This
  controls whether only the averaged rates for the tips, nodes or both
  should be returned. If tips or nodes is chosen the output will be a
  matrix and if both are returned the output will be a list.}
\item{bound.par.matrix}{A numeric matrix with exactly 5 rows and 2
  columns to set the limits for the model parameters. First column is
  the minimum value and second column is the maximum value (both
  inclusive). The rows are the model parameter in this fixed order:
  turnover, net diversification, speciation, extinction fraction, and
  extinction.}
}
\details{
Provides a data frame model-averaged rates for all possible configurations of the model parameters (i.e., turnover, net.div, speciation, extinction, or extinction.fraction), either for all tips or for all nodes.
As with the plotting function, if you give a single hisse.state object,
it uses that, and the rates account for uncertainty in the marginal
probabilities of the states; if you give it a list of them, it will
model-average the results (it assumes the trees are the same) in
addition to accounting for uncertainty in the marginal probabilities at
tips and nodes.
If 'AIC.weights==NULL', then the models in the set do not require the
'$AIC' element to compute the AICw. Otherwise the function will return
an error message.
}
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 
}
\author{Jeremy M. Beaulieu}
\keyword{plotting}
