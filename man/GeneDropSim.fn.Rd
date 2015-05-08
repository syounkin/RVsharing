\name{GeneDropSim.fn}
\alias{GeneDropSim.fn}

\title{Estimation of the probability of sharing of a rare variant by gene dropping in a pedigree}                         
\description{
Estimates the probability that all subjects in a subset of pedigree members share a rare variant given that it occured in any of them by performing a Monte Carlo simulation of the transmission of the genotypes of the variant from the founders down the pedigree.
}
\usage{
GeneDropSim.fn(trio.list, id, dt.vec, fd.indices, carriers=dt.vec, n = 1e3,
 k = 10, nf = 1)
}
\arguments{
  \item{trio.list}{a list of trio objects encoding the pedigree structure. }
  \item{id}{a vector of identifiers of the pedigree members. }
  \item{dt.vec}{ a vector of identifiers of the subset of pedigree members for which to estimate the sharing probability. Must be a subset of the \code{id} vector.}
  \item{fd.indices}{a vector of the indices of the founders of the pedigree.}
  \item{carriers}{optional character vector of subjects carrying the rare variant among all final descendants. If null, all final descendants are assumed to carry the rare variant.}
  \item{n}{minimal number of gene dropping replicates where the rare variant occurs in at least one member of \code{dt.vec}.}
  \item{k}{this number times \code{n} gives the maximal number of gene dropping replicates.}
  \item{nf}{number of founders introducing the rare variant into the pedigree.}
  }
\value{
  Estimate of the probability that all subjects in a subset of pedigree members share a rare variant given that it occured in any of them 
  }
  \details{ The transmission of the RV down the pedigree from the \code{nf} founders introducing it is simulated according to Mendel's laws.  The events that the variant was observed in any of the subjects from \code{dt.vec} and in all of them are then recorded. The simulation continues until the number of replicates where the RV was observed in any of the subjects from \code{dt.vec} reaches \code{n} or the number of replicates reaches \code{k n}. The RV sharing probability is then estimated as the number of replicates where the RV was observed in all subjects from \code{dt.vec} over \code{n} (or the number of replicates where the RV was observed in any of the subjects when \code{k n} replicates are reached).
    }
\references{
Bureau, A., Younkin, S., Parker, M.M., Bailey-Wilson, J.E., Marazita, M.L., Murray, J.C., Mangold, E., Albacha-Hejazi, H., Beaty, T.H. and Ruczinski, I. (2014) Inferring rare disease risk variants based on exact probabilities of sharing by multiple affected relatives. Bioinformatics, 30(15): 2189-96.  
} 
\seealso{
\code{\link{ped2trio}, \link{GeneDropSimExcessSharing.fn}}
}
\examples{
data(ped.list)
plot(ped.list[[54]])
trio.obj = ped2trio(ped.list[[54]])
         
GeneDropSim.fn( trio.list = trio.obj$object, id=ped.list[[54]]$id, dt.vec = c("40","47"),
 fd.indices = trio.obj$fd.indices, n = 5e3)

# Result should be very close to exact value
RVsharing(ped.list[[54]])
}
\author{Samuel G. Younkin <syounkin@stat.wisc.edu> and Alexandre Bureau <alexandre.bureau@msp.ulaval.ca>}
