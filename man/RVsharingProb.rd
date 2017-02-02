\name{RVsharingProb-class}
\docType{class}

% Class:
\alias{class:RVsharingProb}
\alias{RVsharingProb-class}
\alias{RVsharingProb}

\title{RVsharingProb Class}

\description{An object created by \code{RVsharing}}
\arguments{
  \item{pshare}{probability that the \code{carriers} in the pedigree (or all final
descendants in the pedigree if \code{carriers == NULL}) share a rare variant given that a rare variant has been detected in any subject in the union of the \code{carriers} and the final descendants of the pedigree.}
  \item{iancestors}{Character vector of the IDs of branching individuals (intermediate ancestors): subjects who are ancestors to final descendants through two or more of their children and have ancestors above them in the pedigree. The only exception is that one of the top founders is designated as the last branching individual.}
  \item{desfounders}{List of vectors. Each final descendant has a vector in the list containing the distances to the founders above him.}
  \item{id}{Character vector of subject IDs.}
  \item{dad.id}{Character vector of father IDs for the subjects in \code{id}.}
  \item{mom.id}{Character vector of mother IDs for the subjects in \code{id}.}
  \item{carriers}{Character or numeric vector of subjects carrying the rare variant among all final descendants. May be \code{NULL}.}
  
}

\author{Alexandre Bureau <alexandre.bureau@msp.ulaval.ca>}
