\name{RVsharing}
\alias{RVsharing}
\alias{RVsharing,pedigree,missing,missing,character-method}
\alias{RVsharing,pedigree,missing,missing,missing-method}
\alias{RVsharing,character,character,character,character-method}
\alias{RVsharing,character,character,character,missing-method}
\alias{RVsharing,numeric,numeric,numeric,numeric-method}
\alias{RVsharing,numeric,numeric,numeric,missing-method}
\title{Probability of sharing a rare variant among relatives}
\description{
Computing probability that a rare variant is shared by a set of subjects
in a pedigree using equation (1) of Bureau et al. and a generalization to a subset of subjects.
}
\usage{
RVsharing(data, dad.id, mom.id, carriers)
}
\arguments{
  \item{data}{a \code{pedigree} object or character/numeric vector of
    subject IDs.}
  \item{dad.id}{if data is a vector, character or numeric vector of father IDs. Founders' parents should be coded to NA or 0.}
  \item{mom.id}{if data is a vector, character or numeric vector of mother IDs. Founders' parents should be coded to NA or 0.}
  \item{carriers}{optional character or numeric vector of subjects carrying the rare variant. If null, all final descendants are assumed to carry the rare variant.}
}
\value{ An object of class \code{RVsharingProb}
  }
  \details{

    The function \code{RVsharing} computes the probability that all subjects identified as carriers of a rare variant in the vector \code{carriers} (or all final
descendants in the pedigree if \code{carriers = NULL}) share that rare variant AND the final descendants not included in \code{carriers} do not carry it, given that the rare
variant has been detected in any subject in the union of the carriers and the final descendants of the pedigree. A final descendant is defined as a subject without descendant in the pedigree, it it not necessarily in the youngest generation. Multiple mariages can only involve one of the top
founders. Branching individuals (intermediate ancestors) can have only one spouse. Other restrictions remain on the pedigree structure; the function will exit with an error message if such restriction is violated. All final
descendants must share a common ancestor or couple of ancestors, otherwise an erroneous response may be obtained. 
The function recursively processes branching individuals (intermediate ancestors) from the lowest one in the pedigree to one who is a top founder, applying the formulas in Bureau et al. to compute the terms of P[C1 = ... = Ck = 1], and P[C1 = ... = Cn = 0], where n is the number of final descendants and k=<n. If \code{carriers} enumerates a subset of pedigree members, the function will then compute the probability these carriers share the rare variant AND the final descendants not included in \code{carriers} do not carry it based on the above terms. To obtain the probability that a set of pedigree members carry a rare variant given it was seen in any of the set members (ignoring the carrier status of final descendants not in the set), the pedigree must be trimmed of the other final descendants before calling \code{RVsharing}. Important note: the \code{affected} element of the \code{pedigree} object is ignored by \code{RVsharing}. 
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


