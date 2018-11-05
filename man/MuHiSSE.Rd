\name{MuHiSSE}
\alias{MuHiSSE}
\title{Multicharacter Hidden State Speciation and Extinction}
\description{Sets up and executes a MuHiSSE model (Multicharacter Hidden State Speciation
 and Extinction) on a phylogeny and character distribution.}
\usage{
MuHiSSE(phy, data, f=c(1,1,1,1), turnover=c(1,2,3,4), eps=c(1,2,3,4), 
hidden.states=FALSE, trans.rate=NULL, condition.on.survival=TRUE, 
root.type="madfitz", root.p=NULL, sann=FALSE, sann.its=10000, 
bounded.search=TRUE,  max.tol=.Machine$double.eps^.50, starting.vals=NULL,
turnover.upper=10000, eps.upper=3, trans.upper=100, restart.obj=NULL, 
ode.eps=0)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and
  with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a matrix (or dataframe) with three columns. The first column 
  containing the species names and the second and third containing the binary 
  character information. Character "0" is on column 2 and chracter "1" is on 
  column 3. A value of 0 means character absent and a value of 1 character 
  present. The input of data follows a Pagel model. See 'Details'.}
\item{f}{vector of length 4 with the estimated proportion of extant
  species in 00, 01, 10, and 11 that are included in the phylogeny. A value of 
  c(0.50, 0.25, 0.125, 0.125) means that 50 percent of species in combination '00', 25 
  percent in '01' and 12.5 percent in '10' and '11'. By default all species are assumed 
  to be sampled.}
\item{turnover}{a numeric vector of length equal to 4+(number of
  \code{hidden.states} * 4). A MuSSE model has 4 speciation parameters:
  lambda00, lambda01, lambda10, and lambda11. A MuHiSSE model with 
  one hidden states has 8 speciation parameters: lambda0A, s1A, s01A, s0B, s1B, s01B. 
  And so on. The length of the numeric vector needs to match the number of speciation 
  parameters in the model. See 'Details'.}
\item{eps}{a numeric vector of length equal to4+(number of
  \code{hidden.states} * 4). A MuSSE model has 4 extinction parameters:
  mu00, mu01, mu10, and mu11. A MuHiSSE model with one hidden state has 8 extinction 
  parameters: mu00A, mu01A, mu10A, mu11A, mu00B, mu01B, mu10B, and mu11B. And so on. The 
  length of the numeric vector needs to match the number of extinction parameters in the 
  model. See 'Details'.}
\item{hidden.states}{a logical indicating whether the model includes a
  hidden states. The default is \code{FALSE}.} 
\item{trans.rate}{provides the transition rate model. See function
  \code{TransMatMakerMuHiSSE}.}
\item{condition.on.survival}{a logical indicating whether the likelihood
  should be conditioned on the survival of two lineages and the
  speciation event subtending them (Nee et al. 1994). The default is \code{TRUE}.}
\item{root.type}{indicates whether root summarization follow the procedure described by FitzJohn et al. 2009, \dQuote{madfitz} or Herrera-Alsina et al. 2018, \dQuote{herr_als}.}
\item{root.p}{a vector indicating fixed root state probabilities. The
  default is \code{NULL}.}
\item{sann}{a logical indicating whether a two-step optimization
  procedure is to be used. The first includes a simulate annealing
  approach, with the second involving a refinement using
  \code{subplex}. The default is \code{FALSE}.}
\item{sann.its}{a numeric indicating the number of times the simulated
  annealing algorithm should call the objective function.}
\item{bounded.search}{a logical indicating whether or not bounds should
  be enforced during optimization. The default is is \code{TRUE}.}
\item{max.tol}{supplies the relative optimization tolerance to
  \code{subplex}.}
\item{starting.vals}{a vector of starting values for the diversification parameters
  to be used instead of the default settings.}
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
\details{
This function sets up and executes a multiple state HiSSE model. The model allows up to 8
hidden categories (hidden states A-H), and implements a more efficient means of carrying 
out the branch calculation. Specifically, we break up the tree into carry out all 
descendent branch calculations simultaneously, combine the probabilities based on their 
shared ancestry, then repeat for the next set of descendents. In testing, we've found that 
as the number of taxa increases, the calculation becomes much more efficient. In future 
versions, we will likely allow for multicore processing of these calculations to further 
improve speed. We also note that there is vignette tha describes more details for running
this particular function.

As for data file format, \code{MuHiSSE} expects a three column matrix or data frame, with 
the first column containing the species names and the second and third containing the 
binary character information. Note that the order of the data file and the names in the 
\dQuote{phylo} object need not be in the same order; \code{MuHiSSE} deals with this 
internally. Also, the character information must be coded as 0 and 1, otherwise, the
function will misbehave. However, if the state for a species is unknown for either 
character, a user can specify this with a 2, and the state will be considered maximally 
ambiguous for all relevant character combinations. For example, if character 1 is in state
0, but character 2 is provided a 2, then the program provides a probability of 1 for 00 
and a probability of for 01.

As with \code{hisse}, we employ a modified optimization procedure. In other words, rather
than optimizing birth and death separately, \code{MuHisse} optimizes orthogonal 
transformations of these variables: we let tau = birth+death define "net turnover", and 
we let eps = death/birth define the \dQuote{extinction fraction}. This reparameterization 
alleviates problems associated with overfitting when birth and death are highly 
correlated, but both matter in explaining the diversity pattern. 

To setup a model, users input vectors containing values to indicate how
many free parameters are to be estimated for each of the variables in
the model. This is done using the \code{turnover} and
\code{extinct.frac} parameters. One needs to specify a value for each of
the parameters of the model, when two parameters show the same value,
then the parameters are set to be linked during the estimation of the
model. For example, a MuHiSSE model with 1 hidden state and all free
parameters has \code{turnover = 1:8}. The same model with
turnover rates constrained to be the same for all hidden states has
\code{turnover = c(1,2,3,4,1,2,3,4)}. This same format applies to
\code{extinct.frac}.

The \dQuote{trans.rate} input is the transition model and has an entirely different setup 
than speciation and extinction rates. See \code{TransMatMakerMuHiSSE} function for more 
details. 

For user-specified \dQuote{root.p}, you should specify the probability
for each state combination. If you are doing a hidden model, there will be eight state 
combinations: 00A, 01A, 10A, 11A, 00B, 01B, 10B, 11B. So if you wanted to say the root 
had to be in state 00, and since you do not know the hidden state, you would specify 
\dQuote{root.p = c(0.5, 0, 0, 0, 0.5, 0, 0, 0)}. In other words, the root has a 50\% 
chance to be in one of the states to be 00A or 00B.

For the \dQuote{root.type} option, we are currently maintaining the previous default of 
\dQuote{madfitz}. However, it was recently pointed out by Herrera-Alsina et al. (2018)
that at the root, the individual likelihoods for each possible state should be conditioned
prior to averaging the individual likelihoods across states. This can be set doing 
\dQuote{herr_als}. It is unclear to us which is exactly correct, but it does seem that 
both \dQuote{madfitz} and \dQuote{herr_als} behave exactly as they should in the case of 
character-independent diversification (i.e., reduces to likelihood of tree + likelihood of
trait model). We've also tested the behavior and the likelihood differences are very 
subtle and the parameter estimates in simulation are nearly indistinguishable from the
\dQuote{madfitz} conditioning scheme. We provide both options and encourage users to try 
both and let us know conditions in which the result vary dramatically under the two root
implementations. We suspect they do not. 

Also, note, that in the case of \dQuote{root.type=user} and \dQuote{root.type=equal} are
no longer explicit \dQuote{root.type} options. Instead, either \dQuote{madfitz} or 
\dQuote{herr_als} are specified and the \dQuote{root.p} can be set to allow for custom
root options.  
}
\value{
\code{MuHiSSE} returns an object of class \code{muhisse.fit}. This is a list with 
elements:
\item{$loglik}{the maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample-size.}
\item{$solution}{a matrix containing the maximum likelihood estimates of the model 
parameters.}
\item{$index.par}{an index matrix of the parameters being estimated.}
\item{$f}{user-supplied sampling frequencies.}
\item{$hidden.states}{a logical indicating whether hidden states were included in the 
model.}
\item{$condition.on.surivival}{a logical indicating whether the
  likelihood was conditioned on the survival of two lineages and the speciation event 
subtending them.}
\item{$root.type}{indicates the user-specified root prior assumption.}
\item{$root.p}{indicates whether the user-specified fixed root probabilities.}
\item{$phy}{user-supplied tree}
\item{$data}{user-supplied dataset}
\item{$trans.matrix}{the user-supplied transition matrix}
\item{$max.tol}{relative optimization tolerance.}
\item{$starting.vals}{The starting values for the optimization.}
\item{$upper.bounds}{the vector of upper limits to the optimization search.}
\item{$lower.bounds}{the vector of lower limits to the optimization
  search.}
\item{$ode.eps}{The ode.eps value used for the estimation.}
}
\references{
Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models 
of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 

FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating trait-dependent speciation 
and extinction rates from incompletely resolved phylogenies. Syst. Biol. 58:595-611.

Herrera-Alsina L., van Els P., and Etienne R.S. 2018. Detecting the dependence of 
diversification on multiples traits from phylogenetic trees and trait data. Systematic 
Biology In press.

Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary characters effect on 
speciation and extinction. Syst. Biol. 56:701-710.

Nakov, T., Beaulieu, J.M., and Alverson, A.J. 2018. Freshwater diatoms diversify faster 
than marine in both planktonic and benthic habitats. bioRxiv, 
doi: https://doi.org/10.1101/406165.

Nee S., May R.M., and Harvey P.H. 1994. The reconstructed evolutionary process. Philos. 
Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
}

\author{Jeremy M. Beaulieu}
\keyword{models}
