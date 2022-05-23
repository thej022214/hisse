\name{SummarizeMiSSEGreedy}
\alias{SummarizeMiSSEGreedy}
\title{Summarize Results of the Greedy Algorithm}
\description{This provides an overview of all the examined models, including (optionally) a reconstruction of the diversification parameters across the models.}
\usage{
SummarizeMiSSEGreedy(greedy.result, min.weight=0.01, n.cores=1, recon=TRUE) }
\arguments{
\item{greedy.result}{a list of misse.fit objects, returned from MiSSEGreedy or just a list of misse results more generally}
\item{min.weight}{what proportion of the AICc weight a model must have to be included in the model average}
\item{n.cores}{how many cores to use for parallelization}
\item{recon}{Boolean on whether to do the rate reconstructions}
}

\value{
\code{SummarizeMiSSEGreedy} returns a list of containing a data.frame with an overview of models (overview) and an array with rates (rates).
}

\details{
After doing a MiSSEGreedy run, this function provides an overview of the models. The overview object is a data.frame with columns for the model number, its AICc and related measures, the number of free parameters, how long that model took, and whether it is used to reconstruct the diversification parameters. Whether to include a model or not in a reconstruction is up to the user: including very bad models can lead to a reconstruction that is not very good (the model might have very low weight, but if the parameter estimate is still near infinity, for example, it could have a major impact). By default we use models with a weight of at least 0.01.

The rates object is a 3d array: the dimensions are the tips and internal nodes, the parameters being estimated, and the model(s) being used. For example, if you store the SummarizeMiSSEGreedy() output as summarized_results, then summarized_results$rates[,"extinction.fraction","best"] is the extinction fraction estimates for the tips and internal nodes from the best model.
}
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models
of trait-dependent speciation and extinction. Syst. Biol. 65:583-601.

FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating trait-dependent speciation
and extinction rates from incompletely resolved phylogenies. Syst. Biol. 58:595-611.

Herrera-Alsina, L., P. van Els, and R.S. Etienne. 2018. Detecting the dependence of diversification on multiples traits from phylogenetic trees and trait data. Systematic Biology, 68:317-328.

Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary characters effect on
speciation and extinction. Syst. Biol. 56:701-710.

Nee S., May R.M., and Harvey P.H. 1994. The reconstructed evolutionary process. Philos.
Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
}

\author{Brian C. O'Meara}
\keyword{models}
