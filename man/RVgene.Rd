\name{RVgene}
\alias{RVgene}

\title{Probability of sharing of rare variants in a family sample within a gene}                         
\description{
Computing probability of sharing of rare variants in a family sample within a genomic region such as a gene.}
\usage{
RVgene(ped.mat,ped.listfams,sites,fams,pattern.prob.list,nequiv.list,N.list,type="alleles",minor.allele.vec,precomputed.prob=list(0))
}
\arguments{
  \item{ped.mat}{a data.frame or matrix encoding the pedigree information and genotype data in the standard LINKAGE ped format (see PLINK web site [1]). In fact, only the family ID in the first column, the subject ID in the second column, the affection status in the sixth column and the genotype data starting in the seventh column are used (columns 3 to 5 are ignored). Also, family members without genotype data do not need to appear in this matrix. The genotype of each variant can be coded in two ways, each corresponding to a different value of the \code{type} option: a minor allele count on one column, as returned for example by the \code{\link{genotypeToSnpMatrix}} function (\code{type}="count") or the identity of the two alleles on two consecutive columns (\code{type}="alleles").}
  \item{ped.listfams}{a list of \code{pedigree} objects, one object for each pedigree in \code{ped.mat}. }
  \item{sites}{a vector of the column indices of the variant sites to test in \code{ped.mat}. If the argument \code{fams} is provided, the variant sites are tested in each corresponding family in the \code{fams} vector (a variant present in multiple families must then be repeated for every families where it appears).}
  \item{fams}{ an optional character vector of the names of families in \code{ped.mat} and \code{ped.listfams} carrying the variants listed in the corresponding position in \code{sites}. If missing, the names of the families carrying the minor allele at  in are extracted from \code{ped.mat}.}
  \item{pattern.prob.list}{a list of precomputed rare variant sharing probabilities for all possible sharing patterns in the families in \code{ped.mat} and \code{ped.listfams}.}
  \item{nequiv.list}{ an optional vector of the number of configurations of rare variant sharing by the affected subjects corresponding to the same pattern and probability in \code{pattern.prob.list}. Default is a vector of 1s.}
  \item{N.list}{a vector of the number of affected subjects sharing a rare variant in the corresponding pattern in \code{pattern.prob.list}.}
  \item{type}{an optional character string taking value "alleles" or "count". Default is "alleles".}
  \item{minor.allele.vec}{an optional vector of the minor alleles at each site in the \code{sites} vector. It is not needed if \code{type}="count". If it is missing and \code{type}="alleles", the minor allele is assumed to take the value 2.}
  \item{precomputed.prob}{an optional list of vectors precomputed rare variant sharing probabilities for families in \code{ped.mat} and \code{ped.listfams}. The vectors represent probabilities for all the possible values of \code{N.list} for the corresponding family (one probability per value of \code{N.list}).}
  }
\value{
A list with items:
\item{p}{P-value of the exact rare variant sharing test allowing for sharing by a subset of affected subjects.}
\item{pall}{P-value of the exact rare variant sharing test requiring sharing by all affected subjects.}
	  }
  \details{ The function extacts the carriers of the minor allele at each entry in \code{sites} in each family where it is present in \code{ped.mat} (or in the families specified in \code{fams} if that argument is specified). It then computes exact rare variant sharing probabilities in each family for each variant by calling \code{RVsharing}. If multiple rare variants are seen in the same family, the smallest sharing probability among all rare variants is retained. 
  	The joint rare variant sharing probability over all families is obtained as the product of the family-specific probabilities. The p-value of the test allowing for sharing by a subset of affected subjects over the rare variants in the genomic region is then computed as the sum of the probabilities of the possible combinations of sharing patterns among all families with a probability less than or equal to the observed joint probability and a total number of carriers greater than or equal to the sum of the number of carriers in all families, using the values in \code{pattern.prob.list}, \code{nequiv.list} and \code{N.list}.
  	The families where all affected subjects share a rare variant are determined by verifying if the length of the carrier vector equals the maximum value of \code{N.list} for that family. The p-value of the test requiring sharing by all affected subjects is computed by calling \code{get.psubset}.
   }
\references{
[1] http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped	
[2] Bureau, A., Younkin, S., Parker, M.M., Bailey-Wilson, J.E., Marazita, M.L., Murray, J.C., Mangold, E., Albacha-Hejazi, H., Beaty, T.H. and Ruczinski, I. (2014) Inferring rare disease risk variants based on exact probabilities of sharing by multiple affected relatives. Bioinformatics, 30(15): 2189-96.  
} 
\seealso{
\code{\link{RVsharing},\link{get.psubset}}
}
\examples{
data(ped.list)
data(ex.ped.mat)
plot(ped.list[[49]])
# Computation of RV sharing probability for 5 sharing patterns in family 28003
fam28003.pattern.prob = c(RVsharing(ped.list[[49]],carriers=c("36","104","110"))@pshare,
RVsharing(ped.list[[49]],carriers=c("36","104"))@pshare,
RVsharing(ped.list[[49]],carriers=c("104","110"))@pshare,
RVsharing(ped.list[[49]],carriers=c("36"))@pshare,
RVsharing(ped.list[[49]],carriers=c("104"))@pshare)
fam28003.nequiv = c(1,2,1,1,2)
# check that distribution sums to 1
sum(fam28003.pattern.prob*fam28003.nequiv)
fam28003.N = c(3,2,2,1,1)

plot(ped.list[[13]])
# Computation of RV sharing probability for 3 sharing patterns in family 15157
fam15157.pattern.prob = c(RVsharing(ped.list[[13]],carriers=c("402","404","405"))@pshare,
RVsharing(ped.list[[13]],carriers=c("402","404"))@pshare,
RVsharing(ped.list[[13]],carriers=c("402"))@pshare)
fam15157.nequiv = c(1,3,3)
# check that distribution sums to 1
sum(fam15157.pattern.prob*fam15157.nequiv)
fam15157.N = 3:1

# Creating lists
ex.pattern.prob.list = list("15157"=fam15157.pattern.prob,"28003"=fam28003.pattern.prob)
ex.nequiv.list = list("15157"=fam15157.nequiv,"28003"=fam28003.nequiv)
ex.N.list = list("15157"=fam15157.N,"28003"=fam28003.N)
ex.ped.obj = ped.list[c(13,49)]
names(ex.ped.obj) = c("15157","28003")

sites = c(92,119)
minor.allele.vec=c(1,4)

RVgene(ex.ped.mat,ex.ped.obj,sites,pattern.prob.list=ex.pattern.prob.list,nequiv.list=ex.nequiv.list,N.list=ex.N.list,minor.allele.vec=minor.allele.vec)
}
\author{Alexandre Bureau <alexandre.bureau@msp.ulaval.ca>}
