\name{ped2trio}
\alias{ped2trio}

\title{Conversion of a \code{pedigree} object into a list of trio objects}                         
\description{
Creates trio objects representing the pedigree structure contained in the \code{pedigree} object.
}
\usage{
ped2trio(ped)
}
\arguments{
  \item{ped}{a \code{pedigree} object }
  }
\value{
  \item{object}{list of trio objects, one list for each founder couple of the pedigree}
  \item{fd.indices}{Vector of the IDs of the founders of the pedigree}
  }
  \details{

    The function processes the pedigree from the bottom generation to the top, adding new trio objects for each non-founder in a generation, with in their offspring list the trio objects associated to their children. At the top generation, there are as many trio objects as the number of founder couples, and these are returned as a list of trio objects. For now, only founders are allowed to have more than one spouse, non-founders can only have one spouse.
    }
\seealso{
\code{\link{pedigree}}
}    
\examples{
data(ped.list)                                                                                                                                                                                                                                                                                                  
ped2trio(ped.list[[1]])
}
\author{Alexandre Bureau <alexandre.bureau@msp.ulaval.ca>}
