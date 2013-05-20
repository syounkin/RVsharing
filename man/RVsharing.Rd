\name{RVsharing}
\alias{RVsharing}
\title{Main function of package Bureau}
\description{
Computing probability that a rare variant is shared by a set of subjects
in a pedigree using equation 1 of Bureau et al.
}
\usage{
RVsharing(id, dad.id, mom.id)
}
\arguments{
  \item{id}{vector of subject IDs, format unknown}
  \item{dad.id}{vector of father IDs, format unknown}
  \item{mom.id}{vector of mother IDs, format unknown}
}
\value{
  \item{p.share}{no details yet}
  \item{iancestors}{no details yet}
  \item{desfounders}{no details yet}
  \item{id}{no details yet}
  \item{dad.id}{no details yet}
  \item{mom.id}{no details yet}
  }
  \details{

Version 0.2 2013/02/12 Addition of a wrapper for pedigree objects Change
of fdi from a logical vector to a vector of indices, to correct a
problem occuring when an intermediate ancestor is before the final
descendants in the ID list Restriction of currentnonfounders to the
active ones when determining whether there is an intermediate ancestor
Removed currentfounders that are NA when computing how many there are.
At depth 0, Include previous ancestors of only the final descendants
that are active Version 0.3 2013/04/08 Adding the names of the founders
to the desfounders list of degrees between final descendants and
founders Moving the inactivation of final descendants below an
intermediate ancestor and the removal of the spouse of the intermediate
ancestor after the update of the List of distance to founders of each
final descendant to insure inclusion of the spouse in the list Adding
desfounders, id, dad.id and mom.id to the list of returned elements for
further computations with other functions

    The function RV sharing computes the probability that all final
descendants in the pedigree share a rare variant given that a rare
variant has been detected in any one of these final descendants For now,
there can only be one lineage of intermediate ancestors with more than
one child each Multiple mariages can only involve one of the top
founders. Intermediate ancestors can have only one spouse All final
descendants must share a common ancestor or couple of ancestors,
otherwise an erroneous response may be obtained Correction to the
extration of the kth name : names(tmp)[k] }

