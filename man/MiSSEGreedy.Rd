\name{MiSSEGreedy}
\alias{MiSSEGreedy}
\title{Character-free State Speciation and Extinction Searching Greedily}
\description{Executes a set of MiSSE models (Missing State Speciation
 and Extinction) on a phylogeny, varying the number of parameters for turnover and extinction fraction and stopping when models stop being very good.}
\usage{
MiSSEGreedy(phy, f=1, possible.combos = generateMiSSEGreedyCombinations(),
stop.deltaAICc=10, save.file=NULL, n.cores=NULL, chunk.size=NULL,
condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, sann=TRUE,
sann.its=10000, bounded.search=TRUE, max.tol=.Machine$double.eps^.50, starting.vals=NULL,
turnover.upper=10000, eps.upper=3, trans.upper=100, restart.obj=NULL, ode.eps=0)
}
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and
  with internal nodes labeled denoting the ancestral selective regimes.}
\item{f}{the estimated proportion of extant species included in the phylogeny. A value of
  0.50 means that 50 percent of species are contained in the. By default all species are assumed
  to be sampled.}
\item{possible.combos}{data.frame of parameter combinations to try. See 'Details'.}
\item{stop.deltaAICc}{how bad compared to the best does a model have to be to far enough outside we stop looking at even more complex ones?}
\item{save.file}{file to use to save results while the code is running}
\item{n.cores}{how many cores to run this on in parallel}
\item{chunk.size}{how many models to run before checking to make sure there are no improvements. See 'Details'.}
\item{condition.on.survival}{a logical indicating whether the likelihood
  should be conditioned on the survival of two lineages and the
  speciation event subtending them (Nee et al. 1994). The default is \code{TRUE}.}
\item{root.type}{indicates whether root summarization follow the procedure described by FitzJohn et al. 2009, \dQuote{madfitz} or Herrera-Alsina et al. 2018, \dQuote{herr_als}.}
\item{root.p}{a vector indicating fixed root state probabilities. The
  default is \code{NULL}.}
\item{sann}{a logical indicating whether a two-step optimization
  procedure is to be used. The first includes a simulate annealing
  approach, with the second involving a refinement using
  \code{subplex}. The default is \code{TRUE}.}
\item{sann.its}{a numeric indicating the number of times the simulated
  annealing algorithm should call the objective function.}
\item{bounded.search}{a logical indicating whether or not bounds should
  be enforced during optimization. The default is \code{TRUE}.}
\item{max.tol}{supplies the relative optimization tolerance to
  \code{subplex}.}
\item{starting.vals}{a numeric vector of length 3 with starting values for the
  model. Position [1] sets turnover, [2] sets extinction fraction, and
  [3] transition rates between distinct diversification rates.}
\item{turnover.upper}{sets the upper bound for the turnover parameters.}
\item{eps.upper}{sets the upper bound for the eps parameters.}
\item{trans.upper}{sets the upper bound for the transition rate parameters.}
\item{restart.obj}{an object of class  that contains everything to restart an optimization.}
\item{ode.eps}{sets the tolerance for the integration at the end of a
  branch. Essentially if the sum of compD is less than this tolerance,
  then it assumes the results are unstable and discards them. The
  default is set to zero, but in testing a value of 1e-8 can sometimes
  produce stable solutions for both easy and very difficult optimization problems.}
}

\value{
\code{MiSSEGreedy} returns a list of class \code{misse.fit} objects.
}

\details{
See the \code{MiSSE} function for description of the method overall. It requires a set number of hidden state categories, but finding the best number of categories can be hard. For example, one could 1 to 26 different turnover categories and 1 to 26 possible extinction fraction categories. For most cases, we suspect that it makes sense to have the number of extinction fraction categories either equal to the number of turnover categories or set to the same category over the tree, but there are actually a lot of possibilities: have turnover=c(1,2,3) and eps=c(1,2,1), for example, or turnover=c(1,1,1) and eps=c(1,2,3). This uses the generateMiSSEGreedyCombinations function to generate a very large set of these possible models, then runs them in increasing complexity. By default, it stops when the models stop getting being reasonable in AICc. This is NOT where the models stop being significant (if you're looking for significance, note you don't get that from AICc), but where they probably will not contribute much to the model averaged parameter estimates and so may not be worth the bother of searching further. However, there's no guarantee that this is a wise decision: you could stop with 10 turnover rates, and 11 and 12 are far worse for AICc, but it could be that 13 turnover rates are much better than anything else (for that matter, the best number of turnover rates could be 42, even though MiSSe's current code can only go up to 26). And of course, in reality, the truth is that there is an infinite set of hidden rate parameters: a passing cloud blocks a tiny bit of energy for photosynthesis, so for that moment the rate of extinction for plants underneath the shadow is a tiny bit higher. You should not be using this to test hypotheses about how many hidden factors there are, but rather as a way to get a good enough model to estimate rates on a tree.

You can change how quickly the function stops trying new models with stop.deltaAICc. Once it runs a chunk of models where the best gets a model that is at least stop.deltaAICc worse than the *current* best model, it stops running new models. Since this is based on current best AICc, and we start with simple models, there's an asymmetry: a terrible model with no rate variation is always included, but a slightly less terrible model with 26 turnover rates might never be run.

This works MUCH faster when run in parallel. To do so, set n.cores to the number of available cores on your machine (for example, for a modern Mac laptop, this would be 4). An easy way to do this automatically is to use n.cores=parallel::detectCores(). The default is one core; if running this function on a cluster, a default of parallel::detectCores() might take over an entire compute node when you're supposed to be using just one core, and get you in trouble. Since we run many models, the most natural approach is to run one model per core, see if at least one of the models are still ok, then send out the next models out to all the cores. Setting n.cores to the number of parallel jobs you want, and leaving chunk.size set to NULL, will do this. However, this is slightly inefficient -- the odds are that some cores will finish earlier than others, and will be waiting until all finish. So a different approach is to set a chunk.size greater than n.cores -- it will still use no more than n.cores at a time, but once one model in the set finishes it will send off the next until all models in the chunk are run. This keeps the computer even busier, but then it won't stop to check to make sure the models are still feasible as often, and it only saves intermediate results after each chunk of models finishes. Our recommendation is to use n.cores=parallel::detectCores() if you're on a machine where you can use all the cores and leave chunk.size unspecified (so it will also default to n.cores), but it's up to you.

After every chunk of models are done, this function will display the status: what models have been run, what the likelihoods and AICs are, etc. It will also predict how long future runs will check (based on a linear regression between the number of free parameters and log(minutes to run)). These are just estimates based on the runs so far, but it's a stochastic search and can take more or less time.

Saving output while this goes is highly recommended -- you can see how things are going and salvage something if a run fails (then, start again, deleting the models that worked from possible.combos and USING A DIFFERENT FILE TO SAVE TO SO YOU DON'T OVERWRITE THE OLD ONE. You can so this by giving a file name (including path, if you want) as the save.file argument. This will save the list of finished models (misse.list) and the possible.combos data.frame with additional information.
}



\examples{
\donttest{
  library(ape)
  phy <- rcoal(50)
  potential.combos <- generateMiSSEGreedyCombinations(max.param=4, vary.both=TRUE)

  #run on one core:
  result <- MiSSEGreedy(phy, potential.combos, n.cores=1, f=1) 

  #run on all available cores: (commented out because CRAN does not allow this in examples)
  #result <- MiSSEGreedy(phy, potential.combos, n.cores=parallel::detectCores(), f=1)



  AICc <- unlist(lapply(result, "[[", "AICc"))
  deltaAICc <- AICc-min(AICc)
  print(length(result))

  # Yes, [[ is a function you can lapply, and the name of elements within each list
  # object can be arguments. Life is awesome.
}
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
