\name{GeoHiSSE_Plus}
\alias{GeoHiSSE_Plus}
\title{Hidden Geographic State Speciation and Extinction (2+ areas)}
\description{Sets up and executes a GeoHiSSE model (Hidden Geographic State Speciation
 and Extinction) on a phylogeny and character distribution with more
 than 2 endemic areas.}
\usage{
GeoHiSSE_Plus(phy, data, f=c(1,1,1), speciation=c(1,2,3), extirpation=c(1,2), 
hidden.areas=FALSE, trans.rate=NULL, assume.cladogenetic=TRUE, 
condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, sann=FALSE,
sann.its=10000, bounded.search=TRUE,  max.tol=.Machine$double.eps^.50,
mag.san.start=0.5, starting.vals=NULL, speciation.upper=1000, extirpation.upper=1000, 
trans.upper=100, ode.eps=0)
}
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and
  with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a matrix (or dataframe) with 1's and 0's showing presence
  and absence of species on three areas, respectivelly. Each column is
  an area and each row has the distribution information for each
  species. Column names can be used to specify the names of the
  ranges. See examples.}
\item{f}{numeric vector of length 6 with the estimated proportion of extant
  species in each of the areas. Values need to sum to 1. If the areas on
  columns of data are 1, 2, and 3. Then the order of the vector 'f' are
  areas 1, 2, 3, 1+2, 1+3, 2+3.}
\item{speciation}{a numeric vector of length equal to 6+(number of
  \code{hidden.areas} * 6). A GeoSSE model with 3 areas has 6 speciation
  parameters (in this order): s1, s2, s3, s12, s21, and s23. A GeoHiSSE_Plus
  model with 3 areas and one hidden state has 12 speciation parameters:
  s1A, s2A, s3A, s12A, s13A, s23A, s1B, s2B, s3B, s12B, s13B, s23B, and
  so on. The length of the numeric vector needs to match the number of
  speciation parameters in the model. See 'Details'.}
\item{extirpation}{a numeric vector of length equal to 3+(number of
  \code{hidden.areas} * 3). A GeoSSE model with 3 areas has 3
  extirpation parameters: x1, x2, and x3. A GeoHiSSE_Plus model with 3 areas
  one hidden rate has 6 extirpation parameters: x1A, x2A, x3A, x1B, x2B,
  and x3B. And so on. The length of the numeric vector needs to match
  the number of extirpation parameters in the model. See 'Details'.}
\item{hidden.areas}{a logical indicating whether the model includes a
  hidden area. The default is \code{FALSE}.} 
\item{trans.rate}{provides the transition rate model. See function
  \code{TransMatMakerGeoHiSSE}.}
\item{assume.cladogenetic}{assumes that cladogenetic events occur at
  nodes. The default is \code{TRUE}.}
\item{condition.on.survival}{a logical indicating whether the likelihood
  should be conditioned on the survival of two lineages and the
  speciation event subtending them (Nee et al. 1994). The default is \code{TRUE}.}
\item{root.type}{indicates whether root summarization follow the
  procedure described by FitzJohn et al. 2009, \dQuote{madfitz} or
  Herrera-Alsina et al. 2018, \dQuote{herr_als}. }
\item{root.p}{a vector of length 6 indicating fixed root state probabilities. The
  default is \code{NULL}. Order of the areas is s1, s2, s3, s12, s21, and s23.}
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
\item{mag.san.start}{Sets the extinction fraction to estimate the starting values 
  for the diversification parameters. The equation used is based on Magallon and 
  Sanderson (2001), and follows the procedure used in the original GeoSSE 
  implementation.}
\item{starting.vals}{a vector of starting values for the diversification parameters
  to be used instead of the default settings.}
\item{speciation.upper}{sets the upper bound for the speciation parameters.}
\item{extirpation.upper}{sets the upper bound for the extirpation parameters.}
\item{trans.upper}{sets the upper bound for the transition rate parameters.}
\item{ode.eps}{sets the tolerance for the integration at the end of a
  branch. Essentially if the sum of compD is less than this tolerance,
  then it assumes the results are unstable and discards them. The
  default is set to zero, but in testing a value of 1e-8 can sometimes
  produce stable solutions for both easy and very difficult optimization problems.}
}
\details{
This function sets up and executes the GeoHiSSE_Plus model. The model closely
follows \code{diversitree}, although here we employ modified
optimization procedures. As for data file format, \code{GeoHiSSE_Plus}
expects a two column matrix or data frame, with the first column
containing the species names and the second containing the are
information. The area information need to be in the format of three
numbers: 0 for area '01', 1 for area '0', and 2 for '1'. Please note that 
the code for the areas here differ from the
\code{make.geosse} function of package \code{diversitree}.
The order of the data file and the names in the
\dQuote{phylo} object need not be in the same order; \code{hisse} deals
with this internally. Also, the character information MUST be 0,
1, or 2, otherwise, the function will return an error message.

To setup a model, users input vectors containing values to indicate how
many free parameters are to be estimated for each of the variables in
the model. This is done using the \code{speciation} and
\code{extirpation} parameters. One needs to specify a value for each of
the parameters of the model, when two parameters show the same value,
then the parameters are set to be linked during the estimation of the
model. For example, a GeoHiSSE_Plus model with 1 hidden area and all free
parameters has \code{speciation = 1:6}. The same model with
speciation rates constrained to be the same for all hidden areas has
\code{speciation = c(1,2,3,1,2,3)}. This same format applies to
\code{extirpation}. Please note that GeoHiSSE_Plus currently works with up to
4 hidden areas. The most complex model would be \code{speciation = 1:15}
and \code{extirpation = 1:10}.

Once the model is specified, the parameters can be estimated using the
subplex routine (default), or use a two-step process (i.e., sann=TRUE)
that first employs a stochastic simulated annealing procedure, which is
later refined using the subplex routine.

The \dQuote{trans.rate} input is the transition model and has an
entirely different setup than speciation and extirpation rates.
See \code{TransMatMakerGeoHiSSE} function for more details. 

For user-specified \dQuote{root.p}, you should specify the probability
for each area. If you are doing a hidden model, there will be six areas:
0A, 1A, 01A, 0B, 1B, 01B. So if you wanted to say the root had to be in
area 0 (widespread distribution), you would specify \dQuote{root.p =
  c(0.5, 0, 0, 0.5, 0, 0)}. In other words, the root has a 50\% chance
to be in one of the areas 0A or 0B.

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
\code{GeoHiSSE_Plus} returns an object of class \code{geohisse.fit}. This is a list with 
elements:
\item{$loglik}{the maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample-size.}
\item{$solution}{a matrix containing the maximum likelihood estimates of the model parameters.}
\item{$index.par}{an index matrix of the parameters being estimated.}
\item{$f}{user-supplied sampling frequencies.}
\item{$hidden.areas}{a logical indicating whether hidden areas were
  included in the model.}
\item{$assume.cladogenetic}{a logical indicating whether cladogenetic
  events were allowed at nodes.}
\item{$condition.on.surivival}{a logical indicating whether the
  likelihood was conditioned on the survival of two lineages and the
  speciation event subtending them.}
\item{$root.type}{indicates the user-specified root prior assumption.}
\item{$root.p}{indicates whether the user-specified fixed root probabilities.}
\item{$timeslice}{indicates whether the user-specified timeslice that split the tree.}
\item{$phy}{user-supplied tree}
\item{$data}{user-supplied dataset}
\item{$trans.matrix}{the user-supplied transition matrix}
\item{$max.tol}{relative optimization tolerance.}
\item{$starting.vals}{The starting values for the optimization.}
\item{$upper.bounds}{the vector of upper limits to the optimization search.}
\item{$lower.bounds}{the vector of lower limits to the optimization search.}
\item{$ode.eps}{The ode.eps value used for the estimation.}
}
\references{
Caetano, D.S., B.C. O'Meara, and J.M. Beaulieu. 2018. Hidden state models improve state-dependent diversification approaches, including biogeographic models. Evolution, https://doi.org/10.1111/evo.13602

Beaulieu, J.M, and B.C. O'Meara. 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Syst. Biol. 65:583-601. 

FitzJohn R.G., W.P. Maddison, and S.P. Otto. 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Syst. Biol. 58:595-611.

Goldberg, E. E., L. T. Lancaster, and R. H. Ree. 2011. Phylogenetic Inference of Reciprocal Effects between Geographic Range Evolution and Diversification. Syst. Biol. 60:451-465.

Maddison W.P., P.E. Midford, and S.P. Otto. 2007. Estimating a binary characters effect on speciation and extinction. Syst. Biol. 56:701-710.

Nee S., R.M. May, and P.H. Harvey. 1994. The reconstructed evolutionary process. Philos. Trans. R. Soc. Lond. B Biol. Sci. 344:305-311.
}

\author{Jeremy M. Beaulieu}
\keyword{models}
