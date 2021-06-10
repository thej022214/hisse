\name{ProcessSimSample}
\alias{ProcessSimSample}
\title{Process fossilized birth-death sample of lineages}
\description{Takes output from \code{GetFossils} and places extinct lineages in the tree and produces a properly formatted table of fossil samples from lineages that survive to the present (k fossils).}
\usage{
ProcessSimSample(phy, f)
}       
\arguments{
\item{phy}{a complete phylogeny that includes both extant and extinct lineages.} 
\item{f}{the table of sampled lineages under the fossilized birth-death process obtained from \code{GetFossils}.}
}
\value{
Returns an object that contains:
\item{$phy}{a tree of extant taxa and sampled extinct tips.}
\item{$k.samples}{a table of k fossils sampled, formatted to be input into \code{MiSSE}, \code{hisse}, \code{MuHiSSE}.}
}
