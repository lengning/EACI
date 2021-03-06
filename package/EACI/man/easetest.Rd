% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/easetest.R
\name{easetest}
\alias{easetest}
\title{Calculate DAVID EASE score for gene set enrichment analysis}
\usage{
easetest(scores, lib, idtype=c("ENTREZID", "SYMBOL"),library.loc=NULL,sets = c("GO","KEGG"),
locallist = NULL,reduce = NULL,iter=10,minsetsize=10,annotate = TRUE)
}
\arguments{
\item{scores}{numeric vector for gene scores. names should be gene names (IDs)}

\item{lib}{character string, name of data package corresponding to
organism; for example, "org.Hs.eg".}

\item{idtype}{idtype could be either \code{"ENTREZID"} (default) or
\code{"SYMBOL"}. It should match names of the scores vector.}

\item{library.loc}{library location}

\item{sets}{character string, describing the collection of sets. \code{"GO"} (default) or \code{"KEGG"}}

\item{locallist}{list contains in-house gene sets of interest. Default is NULL.
Each element in the list represents a gene set. Each element should be a character
string containing genes' entrez IDs (if \code{idtype="ENTREZID"i}) or
gene symbols (if \code{idtype="SYMBOL"}). Element names will be used as
set names. If \code{locallist} is not NULL, \code{universe="local"} will be disabled.}

\item{reduce}{function used in \code{collapse}.}

\item{iter}{number of iterations when running EM}

\item{minsetsize}{minimal set size to consider; sets with size less than minsetsize will be
ignored during the calculation.}

\item{annotate}{logical, whether to include set names in the output}
}
\value{
sorted list of gene sets; summary statistics
}
\description{
Calculate DAVID EASE score for gene set enrichment analysis
}
\examples{
scores <- rnorm(20)
}
\author{
Ning Leng
}

