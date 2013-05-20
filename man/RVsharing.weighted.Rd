\name{RVsharing.weighted}
\alias{RVsharing}
\title{Main function of package Bureau}
\description{
The function RV sharing computes the probability that all final
descendants in the pedigree share a rare variant given that a rare
variant has been detected in any one of these final descendants
}
\usage{
RVsharing.weighted(id, dad.id, mom.id,relfounders,phi)
}
\arguments{
  \item{id}{vector of subject IDs, format unknown}
  \item{dad.id}{vector of father IDs, format unknown}
  \item{mom.id}{vector of mother IDs, format unknown}
  \item{relfounders}{vector specifying the ids of the founders who are
related, format unkown}
  \item{phi}{either a scalar giving the constant kinhsip coefficient
between the founders in relfounders, or a square symmetric matrix of
kinship coefficients for every pair of founders in relfounders}
}
\details{
For now, there can only be one lineage of intermediate ancestors with
more than one child each Multiple mariages can only involve one of the
top founders.

Intermediate ancestors can have only one spouse

All final descendants must share a common ancestor or couple of
ancestors, otherwise an erroneous response may be obtained

Correction to the extration of the kth name : names(tmp)[k]

The function RVsharing.weighted computes the probability that all
sequenced subjects inherit the variant from a single ancestor (num)
and the probability that no sequenced subject inherited the variant
from any of the founders given that a single copy of the variant is
present in the founders (p0), weighted by the probability that each
founder was the only one to introduce the variant.
}
