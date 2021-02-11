\name{PruneRedundantModels}
\alias{PruneRedundantModels}
\title{Prune Redundant Models}
\description{Removes models that fit equivalently well to models that are as simple or simpler.}
\usage{
PruneRedundantModels(..., precision=1e-5)
}
\arguments{
\item{...}{A list of models (such as from MiSSEGreedy or just putting individual hisse runs in a list)}
\item{precision}{How different the log likelihoods can be to assume they're equivalent models)}
}
\details{
Imagine two models:
Model 1: parameters a and b allowed to vary
Model 2: parameters a set to equal b 

If the maximum likelihood estimate for a and b in model 1 are that they have the same value, it's essentially the same as including model 2 twice, and so other models have lower weight in consequence. This comes up with hidden state models where some hidden states are in extremely low frequencies throughout: they have more parameters than simpler models but the likelihood is essentially the same. Burnham and Anderson (2003, page 342ff) recommend removing redundant models in cases like this where two models are effectively identical. 

This function does this by ordering models by the number of free parameters, then deleting models that have log likelihoods nearly exactly equal (to the precision set above) to a model that is as simple or simpler. This is not a perfect solution -- it's possible that models could have very different parameters and have the same likelihood, for example. Also the "as simple" means that if there are two models with equal numbers of parameters, whichever is first gets saved and the second is removed. 
}
\author{Brian O'Meara}
