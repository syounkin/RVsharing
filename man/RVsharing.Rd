\name{RVsharing}
\alias{RVsharing}
\alias{RVsharing,pedigree,missing,missing,missing-method}
\alias{RVsharing,missing,character,character,character-method}
\alias{RVsharing,missing,numeric,numeric,numeric-method}
\title{Probability of sharing a rare variant among relatives}
\description{
Computing probability that a rare variant is shared by a set of subjects
in a pedigree using equation (1) of Bureau et al.
}
\usage{
RVsharing(ped)
## Default S3 method:
RVsharing(id, dad.id, mom.id)
}
\arguments{
  \item{ped}{a \code{pedigree} object}
  \item{id}{character or numeric vector of subject IDs}
  \item{dad.id}{character or numeric vector of father IDs. Founders' parents should be coded to NA or 0}
  \item{mom.id}{character or numeric vector of mother IDs. Founders' parents should be coded to NA or 0}
}
\value{
  \item{p.share}{probability that all final
descendants in the pedigree share a rare variant given that a rare
variant has been detected in any one of these final descendants}
  \item{iancestors}{Vector of the IDs of branching individuals (intermediate ancestors): subjects who are ancestors to final descendants through two or more of their children and have ancestors above them in the pedigree. The only exception is that one of the top founders is designated as the last branching individual.}
  \item{desfounders}{List of vectors. Each final descendant has a vector in the list containing the distances to the founders above him.}
  \item{id}{vector of subject IDs passed as argument}
  \item{dad.id}{vector of father IDs passed as argument}
  \item{mom.id}{vector of mother IDs passed as argument}
  }
  \details{

    The function RV sharing computes the probability that all final
descendants in the pedigree share a rare variant given that a rare
variant has been detected in any one of these final descendants. For now,
there can only be one lineage of branching individuals (intermediate ancestors) with more than
one child each. Multiple mariages can only involve one of the top
founders. Branching individuals can have only one spouse. All final
descendants must share a common ancestor or couple of ancestors,
otherwise an erroneous response may be obtained. 

The function recursively processes branching individuals (intermediate ancestors) from the lowest one in the pedigree to the one who is a top founder, applying the formulas in Bureau et al. to compute the terms of P[C1 = ... = Cn = 1] and P[C1 = ... = Cn = 0] involving each branching individual.
}
\examples{
plot(ped.list[[1]])
RVsharing(ped.list[[1]])
}
\author{Alexandre Bureau <alexandre.bureau@msp.ulaval.ca>}
\references{Bureau, A et al. Sharing of rare variants by affected relatives: building evidence for causal variants based on exact sharing probabilities. Submitted.
}

