\name{GeneDrop-method}
\docType{methods}

\alias{GeneDrop}
\alias{GeneDrop,Trio,numeric-method}

\title{Gene Dropping method}

\description{Simulate transmission of alleles from the parents to their offsrping according to Mendel's laws.}

\usage{
\S4method{GeneDrop}{Trio,numeric}(trio,geno.vec)
}

\arguments{
  \item{trio}{
    A Trio object. 
  }
  \item{geno.vec}{
    A named genotype vector.
  }
}
\details{Simulate transmission of alleles from the parents to their offsrping in a \code{trio} object according to Mendel's laws.}
\value{A numeric genotype vector.}
\author{Samuel G. Younkin <syounkin@jhsph.edu>}
