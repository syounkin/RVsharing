\name{ComputeKinshipPropCoef}
\alias{ComputeKinshipPropCoef}

\title{Ratio of excess kinship among descendants over mean kinship among founders}
\description{
Computes, for each pair of final descendants in the pedigree structure contained in the \code{pedigree} object, the ratio of the difference between the inferred and expected kinship coefficient for the pair over the mean kinship among founders. 
}
\usage{
ped2trio(ped)
}
\arguments{
  \item{ped}{a \code{pedigree} object.}
  }
\value{
A symmetric matrix of ratios for all pair of final descendants in the pedigree structure contained in the \code{pedigree} object.
}
\detail{ 
The ratio for each pair of final descendants is computed using equation (A1) of Bureau et al. Dividing the difference between the inferred and expected kinship coefficient for each pair by this ratio gives a pair-specific estimate of the mean kinship among founders, which can then be averaged over all pairs of final descendants from the same population to obtain a global estimate of the mean kinship among founders.
}
\references{
Bureau, A., Younkin, S., Parker, M.M., Bailey-Wilson, J.E., Marazita, M.L., Murray, J.C., Mangold, E., Albacha-Hejazi, H., Beaty, T.H. and Ruczinski, I. (under review) Inferring rare disease risk variants based on exact probabilities of sharing by multiple affected relatives.  
} 
\seealso{
\code{\link{pedigree}}
}    
\examples{
data(ped.list)
plot(ped.list[[1]])
ComputeKinshipPropCoef(ped.list[[1]])
}
\author{Alexandre Bureau <alexandre.bureau@msp.ulaval.ca>}
}
