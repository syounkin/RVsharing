\name{get.psubset}
\alias{get.psubset}

\title{Probability of sharing of rare variants in a subset of families}                         
\description{
Computing probability of sharing of rare variants in a subset of families where rare variants are seen based on precomputed family-specific rare variant sharing probabilities.}
\usage{
get.psubset(vec,not,pshare.data)
}
\arguments{
  \item{vec}{a vector of names of all families where a variant is seen.}
  \item{not}{a vector of names of families where not all affected subjects share the rare variant (RV). }
  \item{pshare.data}{a data frame with at least the two following columns:
	    \code{pshare} : vector of RV sharing probabilities;
	    \code{ped.tocompute.vec} : vector of names of the families whose sharing probability is contained in pshare. 
							 The names in the arguments \code{vec} and \code{not} must be found in \code{ped.tocompute.vec}.}
		}
\value{
P-value of the exact rare variant sharing test requiring sharing by all affected subjects.
	  }
  \details{  All the subsets of families of size equal or inferior to the lenght of \code{not} are created, and the joint  probability of each such subset not sharing a rare variant and the remaining families sharing a rare variant is obtained as the product of the family-specific rare variant sharing probabilities or its complement. The function then sums the pattern probabilities inferior or equal to the probability of the observed pattern of the \code{not} families not sharing a rare variant and the remaining families sharing a rare variant.
 	  	   }
\references{
Bureau, A., Younkin, S., Parker, M.M., Bailey-Wilson, J.E., Marazita, M.L., Murray, J.C., Mangold, E., Albacha-Hejazi, H., Beaty, T.H. and Ruczinski, I. (2014) Inferring rare disease risk variants based on exact probabilities of sharing by multiple affected relatives. Bioinformatics, 30(15): 2189-96.  
} 
\seealso{
\code{\link{RVsharing},\link{RVgene}}
}
\examples{
data(ped.list)
# Setting up pshare.data data frame
pshare.data=data.frame(pshare=sapply(ped.list[1:39],function(ped) RVsharing(ped)@pshare),ped.tocompute.vec=sapply(ped.list[1:39],function(ped) ped$famid[1]))
vec = c("15008","15053","15157")
get.psubset(vec,not="15008",pshare.data)
get.psubset(vec,not=c("15008","15053"),pshare.data)
}
\author{Alexandre Bureau <alexandre.bureau@msp.ulaval.ca>}
