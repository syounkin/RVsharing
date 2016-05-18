\name{ComputeKinshipPropCoef-method}
\docType{methods}
\alias{ComputeKinshipPropCoef}
\alias{ComputeKinshipPropCoef,RVsharingProb-method}

\title{Ratio of excess kinship among descendants over mean kinship among founders}
\description{
Computes, for each pair of final descendants in the pedigree structure contained in the \code{pedigree} object, the ratio of the difference between the inferred and expected kinship coefficient for the pair over the mean kinship among founders. 
}
\usage{
\S4method{ComputeKinshipPropCoef}{RVsharingProb}(obj)
}
\arguments{
  \item{obj}{obj is a \code{RVsharingProb} object returned by the function \code{RVsharing}}
  }
\value{
A symmetric matrix of ratios for all pair of final descendants in the pedigree structure contained in the \code{pedigree} object.
}
\details{ 
The ratio for each pair of final descendants is computed using equation (A1) of Bureau et al. Dividing the difference between the inferred and expected kinship coefficient for each pair by this ratio gives a pair-specific estimate of the mean kinship among founders, which can then be averaged over all pairs of final descendants from the same population to obtain a global estimate of the mean kinship among founders.
}
\references{
Bureau, A., Younkin, S., Parker, M.M., Bailey-Wilson, J.E., Marazita, M.L., Murray, J.C., Mangold, E., Albacha-Hejazi, H., Beaty, T.H. and Ruczinski, I. (2014) Inferring rare disease risk variants based on exact probabilities of sharing by multiple affected relatives. Bioinformatics, 30(15): 2189-96, doi:10.1093/bioinformatics/btu198.  
} 
\seealso{
\code{\link{pedigree}}
}    
\examples{
data(ped.list)
plot(ped.list[[1]])
obj = RVsharing(ped.list[[1]])
ComputeKinshipPropCoef(obj)
}
\author{Alexandre Bureau <alexandre.bureau@msp.ulaval.ca>}

