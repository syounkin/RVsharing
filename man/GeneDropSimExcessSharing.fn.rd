\name{GeneDropSimExcessSharing.fn}
\alias{GeneDropSimExcessSharing.fn}

\title{Estimation of the probability of sharing of a rare variant by gene dropping in a pedigree}                         
\description{
Estimates the probability that all subjects in a subset of pedigree members share a rare variant (RV) given that it occured in any of them by performing a Monte Carlo simulation of the transmission of the genotypes of the variant from the founders down the pedigree.
}
\usage{
GeneDropSimExcessSharing.fn(trio.list, id, dt.vec, fd.indices, phihat, RVfreq, ord=5, n = 1e3, k = 10, nf = 1)
}
\arguments{
  \item{trio.list}{a list of trio objects encoding the pedigree structure. }
  \item{id}{a vector of identifiers of the pedigree members. }
  \item{dt.vec}{ a vector of identifiers of the subset of pedigree members for which to estimate the sharing probability. Must be a subset of the \code{id} vector.}
  \item{fd.indices}{a vector of the indices of the founders of the pedigree.}
  \item{phihat}{a vector of values of the mean kinship coefficient between founders. Must be non-negative.} 
  \item{RVfreq}{frequency of the variant in the population (optional). When missing, the variant frequency tends to 0.}
  \item{ord}{order of the polynomial approximation of the number of distinct alleles among pedigree founders.}
  \item{n}{minimal number of gene dropping replicates where the rare variant occurs in at least one member of \code{dt.vec}.}
  \item{k}{this number times \code{n} gives the maximal number of gene dropping replicates.}
  \item{nf}{number of founders introducing the rare variant into the pedigree.}
  }
\value{
  Estimate of the probability that all subjects in a subset of pedigree members share a rare variant given that it occured in any of them 
  }
  \details{ The function performs the following steps. It first determines the probability w that the RV was introduced only once in the pedigree and its complement 1-w that it was introduced twice based on the mean kinship among founders \code{phi.hat}. It then samples an indicator variable of whether one or two copies of the RV were introduced into the family with probability w and 1-w respectively. In practice, this is done by sampling the number of distinct alleles a from an approximate distribution derived from \code{phi.hat}, then sampling the RV among the $a$ alleles. The RV is introduced twice if it is one of the first 2n_f - a alleles, and introduced once otherwise. If it is introduced twice, the pair of founders introducing the RV is sampled with equal probability for all pairs. If it is introduced once, the sole founder introducing it is sampled instead. Then the transmission of the RV down the pedigree from the one or two founders introducing it is simulated according to Mendel's laws.
   The events that the variant was observed in any of the subjects from \code{dt.vec} and in all of them are then recorded. The simulation continues until the number of replicates where the RV was observed in any of the subjects from \code{dt.vec} reaches \code{n} or the number of replicates reaches \code{k n}. The RV sharing probability is then estimated as the number of replicates where the RV was observed in all subjects from \code{dt.vec} over \code{n} (or the number of replicates where the RV was observed in any of the subjects when \code{k n} replicates are reached).
  }
\references{
Bureau, A., Younkin, S., Parker, M.M., Bailey-Wilson, J.E., Marazita, M.L., Murray, J.C., Mangold, E., Albacha-Hejazi, H., Beaty, T.H. and Ruczinski, I. (under review) Inferring rare disease risk variants based on exact probabilities of sharing by multiple affected relatives.  
} 
\seealso{
\code{\link{ped2trio}, \link{GeneDropSim.fn}}
}
\examples{
data(ped.list)
plot(ped.list[[54]])
trio.obj = ped2trio(ped.list[[54]])
         
GeneDropSimExcessSharing.fn( trio.list = trio.obj$object, id=ped.list[[54]]$id, dt.vec = c("40","47"), fd.indices = trio.obj$fd.indices, phihat=0.005,RVfreq=0.01, ord=5, n = 1e4)

# Result should be higher than exact value under the assumption of a variant frequency tending to 0 and no unknown relationship among founders
RVsharing(ped.list[[54]])$pshare
}
\author{Samuel G. Younkin <syounkin@jhsph.edu> and Alexandre Bureau <alexandre.bureau@msp.ulaval.ca>}
