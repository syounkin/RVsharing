\name{RVsharing}
\alias{RVsharing}
\alias{RVsharing,pedigree,missing,missing-method}
\alias{RVsharing,character,character,character-method}
\alias{RVsharing,numeric,numeric,numeric-method}
\title{Probability of sharing a rare variant among relatives}
\description{
Computing probability that a rare variant is shared by a set of subjects
in a pedigree using equation (1) of Bureau et al.
}
\usage{
RVsharing(data, dad.id, mom.id, carriers=NULL)
}
\arguments{
  \item{data}{a \code{pedigree} object or character/numeric vector of
    subject IDs.}
  \item{dad.id}{if data is a vector, character or numeric vector of father IDs. Founders' parents should be coded to NA or 0.}
  \item{mom.id}{if data is a vector, character or numeric vector of mother IDs. Founders' parents should be coded to NA or 0.}
  \item{carriers}{optional character or numeric vector of subjects carrying the rare variant among all final descendants. If null, all final descendants are assumed to carry the rare variant.}
}
\value{ An object of class \code{RVsharingProb}
  }
  \details{

    The function \code{RVsharing} computes the probability that all subjects identified as carriers of a rare variant in the vector \code{carriers} (or all final
descendants in the pedigree if \code{carriers = NULL}) share that rare variant given that the rare
variant has been detected in any one of the final descendants of the pedigree. Multiple mariages can only involve one of the top
founders. Branching individuals (intermediate ancestors) can have only one spouse. Other restrictions remain on the pedigree structure; the function will exit with an error message if such restriction is violated. All final
descendants must share a common ancestor or couple of ancestors, otherwise an erroneous response may be obtained. 
The function recursively processes branching individuals (intermediate ancestors) from the lowest one in the pedigree to the one who is a top founder, applying the formulas in Bureau et al. to compute the terms of P[C1 = ... = Cn = 1] and P[C1 = ... = Cn = 0], where n is the number of final descendants. If \code{carriers} enumerates a subset of final descendants, the function will then compute the probability these carriers share the rare variant AND the other final descendants do not carry it. To obtain the probability that a subset of final descendants carry a rare variant ignoring the carrier status of other final descendants, the pedigree must be trimmed to remove these other final descendants. Important note: the \code{affected} element of the \code{pedigree} object is ignored by \code{RVsharing} 
}
\references{
Bureau, A., Younkin, S., Parker, M.M., Bailey-Wilson, J.E., Marazita, M.L., Murray, J.C., Mangold, E., Albacha-Hejazi, H., Beaty, T.H. and Ruczinski, I. (2014) Inferring rare disease risk variants based on exact probabilities of sharing by multiple affected relatives. Bioinformatics, 30(15): 2189-96. 
} 
\seealso{
\code{\link{pedigree}, \link{RVsharingProb}, \link{GeneDropSim.fn}, \link{GeneDropSimExcessSharing.fn}}
}    
\examples{
data(ped.list)
plot(ped.list[[1]])
RVsharing(ped.list[[1]])
}
\author{Alexandre Bureau <alexandre.bureau@msp.ulaval.ca>}


