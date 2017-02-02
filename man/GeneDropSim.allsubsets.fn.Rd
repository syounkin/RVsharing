\name{GeneDropSim.allsubsets.fn}
\alias{GeneDropSim.allsubsets.fn}

\title{Estimation of the probability of sharing of a rare variant by gene dropping in a pedigree}                         
\description{
Estimates the probabilities of sharing a rare variant for all subsets of a specified set of pedigree members given that the rare variant occured in any subject in the set by performing a Monte Carlo simulation of the transmission of the genotypes of the variant from the founders down the pedigree.
}
\usage{
GeneDropSim.allsubsets.fn(trio.list, id, dt.vec, fd.indices, carriers=dt.vec, n = 1e3,
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
  Vector of estimates of the probabilities of sharing of a rare variant for all subsets of subjects of a specified set of pedigree members given the rare variant occured in any of them. The names of the vector elements are the names of the subjects in the corresponding subset, separated by semi-columns.
  }
  \details{ The transmission of the RV down the pedigree from the \code{nf} founders introducing it is simulated according to Mendel's laws.  The events that the variant was observed in any of the subjects from \code{dt.vec} and in all subsets of subjects included in \code{carriers} are then recorded. The simulation continues until the number of replicates where the RV was observed in any of the subjects from \code{dt.vec} reaches \code{n} or the number of replicates reaches \code{k n}. The RV sharing probabilities are then estimated as the number of replicates where the RV was observed in each subset over \code{n} (or the number of replicates where the RV was observed in any of the subjects when \code{k n} replicates are reached).
    }
\references{
Bureau, A., Younkin, S., Parker, M.M., Bailey-Wilson, J.E., Marazita, M.L., Murray, J.C., Mangold, E., Albacha-Hejazi, H., Beaty, T.H. and Ruczinski, I. (2014) Inferring rare disease risk variants based on exact probabilities of sharing by multiple affected relatives. Bioinformatics, 30(15): 2189-96, doi:10.1093/bioinformatics/btu198.  
} 
\seealso{
\code{\link{ped2trio}, \link{GeneDropSim.fn}}
}
\examples{
data(ped.list)
plot(ped.list[[54]])
trio.obj = ped2trio(ped.list[[54]])
         
GeneDropSim.allsubsets.fn( trio.list = trio.obj$object, id=ped.list[[54]]$id, dt.vec = c("40","47"),
 fd.indices = trio.obj$fd.indices, n = 5e3)

}
\author{Samuel G. Younkin <syounkin@stat.wisc.edu> and Alexandre Bureau <alexandre.bureau@msp.ulaval.ca>}
