\name{RVgene_allshare}
\alias{RVgene_allshare}

\title{Probability of sharing of rare variants among all affected relatives in a family sample within a gene}                         
\description{
Computing probability of sharing of rare variants among all affected relatives in a family sample within a genomic region such as a gene.}
\usage{
RVgene_allshare(ped.mat,ped.listfams,sites,fams,pshare.vec,
type="alleles",minor.allele.vec,precomputed.prob=list(0),maxdim = 1e9)
}
\arguments{
  \item{ped.mat}{a data.frame or matrix encoding the pedigree information and genotype data in the standard LINKAGE ped format (see PLINK web site [1]). In fact, only the family ID in the first column, the subject ID in the second column, the affection status in the sixth column and the genotype data starting in the seventh column are used (columns 3 to 5 are ignored). Also, family members without genotype data do not need to appear in this matrix. The genotype of each variant can be coded in two ways, each corresponding to a different value of the \code{type} option: a minor allele count on one column, as returned for example by the \code{genotypeToSnpMatrix} function, with missing values coded \code{NA} (\code{type}="count") or the identity of the two alleles on two consecutive columns, with missing values coded 0 (\code{type}="alleles").}
  \item{ped.listfams}{a list of \code{pedigree} objects, one object for each pedigree in \code{ped.mat}. }
  \item{sites}{a vector of the column indices of the variant sites to test in \code{ped.mat}. If the argument \code{fams} is provided, the variant sites are tested in each corresponding family in the \code{fams} vector (a variant present in multiple families must then be repeated for every families where it appears).}
  \item{fams}{ an optional character vector of the names of families in \code{ped.mat} and \code{ped.listfams} carrying the variants listed in the corresponding position in \code{sites}. If missing, the names of the families carrying the minor allele at each position in \code{sites} are extracted from \code{ped.mat}.}
  \item{pshare.vec}{a vector of the probabilities that all affected relatives share a rare variant, for every family. This vector must be named with the family names.}
  \item{type}{an optional character string taking value "alleles" or "count". Default is "alleles".}
  \item{minor.allele.vec}{an optional vector of the minor alleles at each site in the \code{sites} vector. It is not needed if \code{type}="count". If it is missing and \code{type}="alleles", the minor allele is assumed to take the value 2.}
  \item{precomputed.prob}{an optional list of vectors precomputed rare variant sharing probabilities for families in \code{ped.mat} and \code{ped.listfams}. The vectors represent probabilities for all the possible values of \code{N.list} for the corresponding family (one probability per value of \code{N.list}).}
  \item{maxdim}{upper bound on the dimension of the array containing the joint distribution of the sharing patterns for all families in \code{fams} (to avoid running out of memory).}
  }
\value{
P-value of the exact rare variant sharing test requiring sharing by all affected subjects.
	  }
  \details{ The function extracts the carriers of the minor allele at each entry in \code{sites} in each family where it is present in \code{ped.mat} (or in the families specified in \code{fams} if that argument is specified). It then computes exact rare variant sharing probabilities in each family for each variant by calling \code{RVsharing}. If multiple rare variants are seen in the same family, the smallest sharing probability among all rare variants is retained. 
  	The families where all affected subjects share a rare variant are determined by verifying if the length of the carrier vector equals the number of affected subjects that family. The joint rare variant sharing probability over all families is obtained as the product of the family-specific probabilities. The p-value of the test requiring sharing by all affected subjects is computed by calling \code{get.psubset}.
   }
\references{
[1] http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped	
[2] Bureau, A., Younkin, S., Parker, M.M., Bailey-Wilson, J.E., Marazita, M.L., Murray, J.C., Mangold, E., Albacha-Hejazi, H., Beaty, T.H. and Ruczinski, I. (2014) Inferring rare disease risk variants based on exact probabilities of sharing by multiple affected relatives. Bioinformatics, 30(15): 2189-96, doi:10.1093/bioinformatics/btu198. 
} 
\seealso{
\code{\link{RVsharing},\link{get.psubset}}
}
\examples{
data(ped.list)
data(ex.ped.mat)
plot(ped.list[[49]])
plot(ped.list[[13]])

ex.ped.obj = ped.list[c(13,49)]
names(ex.ped.obj) = c("15157","28003")
# Computation of RV sharing probability by all affected relatives in families 28003 and 15157
pshare.vec = c(RVsharing(ped.list[[13]],)@pshare,RVsharing(ped.list[[49]])@pshare)
names(pshare.vec) = c("15157","28003")

sites = c(92,119)
minor.allele.vec=c(1,4)

RVgene_allshare(ex.ped.mat,ex.ped.obj,sites,pshare.vec=pshare.vec,minor.allele.vec=minor.allele.vec)
}
\author{Alexandre Bureau <alexandre.bureau@msp.ulaval.ca>}
