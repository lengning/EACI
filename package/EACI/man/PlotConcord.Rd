% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/PlotConcord.R
\name{PlotConcord}
\alias{PlotConcord}
\title{Generate concordance plot of enrichment results}
\usage{
PlotConcord(Out, N=100, stat="set.mean",minsize=10)
}
\arguments{
\item{Out}{output of allez(), easefun() or eacifun()}

\item{N}{top n gene sets to show}

\item{stat}{sort by which summary statistics}

\item{minsize}{the minimal set size to consider}
}
\value{
a heatmap showing overlapping rates of top N gene sets.
}
\description{
Generate concordance plot of enrichment results
}
\examples{
scores <- rnorm(20)
}
\author{
Ning Leng
}

