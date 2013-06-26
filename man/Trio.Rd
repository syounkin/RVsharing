\name{Trio-class}
\docType{class}

% Class:
\alias{class:Trio}
\alias{Trio-class}
\alias{Trio}

\title{Trio Class}

\description{Foo!}

\author{Samuel G. Younkin <syounkin@jhsph.edu>}

%% % Constructor
%% \alias{Trio,SummarizedExperiment,PedClass-method}

%% % Accessors:
%% \alias{pedigree,Trio-method}
%% \alias{pedigree}
%% \alias{MAF,Trio-method}
%% \alias{MAF}
%% \alias{completeTrios,Trio-method}
%% \alias{completeTrios}
%% \alias{parents,Trio-method}
%% %\alias{parents}
%% \alias{baf}
%% \alias{baf,Trio-method}
%% \alias{cnv}
%% \alias{cnv,Trio-method}
%% \alias{geno}
%% \alias{geno,Trio-method}
%% \alias{logR}
%% \alias{logR,Trio-method}
%% %\alias{ids}
%% %\alias{ids,

%% % Methods
%% %\alias{[,Trio-method]}
%% \alias{[,Trio,ANY,ANY,ANY-method}
%% %\alias{[,Trio,ANY,ANY,ANY]}
%% %\alias{[,Trio,ANY,ANY-method]}
%% %\alias{[,Trio,ANY,ANY]}
%% %\alias{[,Trio-method]}


%% \usage{

%% ## Constructors

%% Trio(se, pedigree)
%% \S4method{Trio}{SummarizedExperiment}(se, pedigree)
%% %\S4method{SummarizedExperiment}{missing}(assays, ...)
%% %\S4method{SummarizedExperiment}{list}(assays, ...)
%% %\S4method{SummarizedExperiment}{matrix}(assays, ...)

%% }
%%   \section{Accessors}{
%%   In the following code snippets, \code{fe} is a Trio
%%   object.
%%   \describe{
%%     \item{}{
%%       \code{pedigree(fe)}: Returns the pedigree object
%%     contained in \code{fe}.  The returned object will have class
%%     PedClass.
%%   }
%%   \item{}{
%%     \code{MAF(fe)}: Returns a numeric vector of minor allele
%%     frequencies.  The frequencies are computed after subsetting the data
%%   to be parents only.  This clearly is problematice if your pedigree
%%   contains multiple generations.  This method is essentially a wrapper
%%   for the The minor allele frequency method of SnpStats.
%% }
%% \item{}{
%% \code{completeTrios(fe)}: Returns an \code{n} by 3 (character or factor)?
%%   matrix with rows that correspond to the parent-offspring trios
%%   contained in \code{fe} such that each of the trio members appears in
%%   the assay data.  \code{n} is the number of ``complete trios.''
%%   }
%% \item{}{
%% \code{parents(fe)}: Returns the parents IDs.  Note that this is
%%   intended to be used with parent-offspring trios.  It will return any
%%   subject that has both parents in the pedigree data.
%%   }
%% \item{}{
%% \code{baf(fe)}: Retrieve matrix of B allele frequencies.  Not yet implemented. 
%%   }
%% \item{}{
%% \code{logr(fe)}: Retrieve matrix of log R ratio values.  Not yet implemented. 
%%   }
%% \item{}{
%% \code{geno(fe)}: Retrieve a matrix of genotypes.
%%   }
%% \item{}{
%% \code{cnv(fe)}:
%%   }
%% \item{}{
%% \code{fe[i,j]}: subset operator
%%   }
%% }
%% }
