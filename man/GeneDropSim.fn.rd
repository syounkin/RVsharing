\name{GeneDropSim.fn}
\alias{GeneDropSim.fn}

\title{Estimation of the probability of sharing of a rare variant by gene dropping in a pedigree}                         
\description{
Estimates the probability that all subjects in a subset of pedigree members share a rare variant given that it occured in any of them by performing a Monte Carlo simulation of the transmission of the genotypes of the variant from the founders down the pedigree.
}
\usage{
GeneDropSim.fn <- function(trio.list, geno.vec, dt.vec, fd.indices, n = 1e3, k = 10, nf = 1)
}
\arguments{
  \item{trio.list}{a list of trio objects encoding the pedigree structure }
  \item{id}{a vector of identifiers of the pedigree members }
  \item{dt.vec}{ a vector of identifiers of the subset of pedigree members for which to estimate the sharing probability. Must be a subset of the \code{id} vector.}
  \item{fd.indices}{a vector of the indices of the founders of the pedigree.}
  \item{n}{minimal number of gene dropping replicates where the rare variant occurs in at least one member of \code{dt.vec}}
  \item{k}{this number times \code{n} gives the maximal number of gene dropping replicates.}
  \item{nf}{number of founders introducing the rare variant into the pedigree}
  }
\value{
  Estimate of the probability that all subjects in a subset of pedigree members share a rare variant given that it occured in any of them 
  }
  \details{ Foo
    }
\examples{
data(ped.list)
plot(ped.list[[54]])
trio.obj = ped2trio(ped.list[[54]])                                                                                                                                                                                                                                                                                                  
GeneDropSim.fn( trio = trio.obj$object, id=ped.list[[54]]$id, dt.vec = c("40","47"), fd.indices = trio.obj$fd.indices, n = 1e4)
# Result should be very close to exact value
RVsharing(ped.list[[54]])$pshare
}
\author{Samuel G. Younkin <syounkin@jhsph.edu> and Alexandre Bureau <alexandre.bureau@msp.ulaval.ca>}
}

