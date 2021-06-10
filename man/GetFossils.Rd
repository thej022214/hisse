\name{GetFossils}
\alias{GetFossils}
\title{Get fossilized birth-death samples}
\description{Samples a simulated tree of extinct and extant lineages according to the preservation rate of the fossilized birth-death process.}
\usage{
GetFossils(phy, psi=0.1)
}       
\arguments{
\item{phy}{a complete phylogeny that includes both extant and extinct lineages.} 
\item{psi}{the preservation rate for which lineages will be sampled under a fossilized birth-death process.}
}
\value{
Returns a table of sampled lineages.
}
