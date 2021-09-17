\name{ProcessSimStrat}
\alias{ProcessSimStrat}
\title{Process fossilized birth-death sample of lineages into stratigraphic intervals.}
\description{Takes output from \code{GetFossils} and places extinct lineages in the tree and produces a properly formatted table of stratigraphic intervals.}
\usage{
ProcessSimStrat(phy, f)
}       
\arguments{
\item{phy}{a complete phylogeny that includes both extant and extinct lineages.} 
\item{f}{the table of sampled lineages under the fossilized birth-death process obtained from \code{GetFossils}.}
}
\value{
Returns an object that contains:
\item{$phy}{a tree of extant taxa and sampled extinct tips.}
\item{$strat.intervals}{a table of stratigraphic intervals, formatted to be input into \code{MiSSE}, \code{hisse}, \code{MuHiSSE}.}
}
