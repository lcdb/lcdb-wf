% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{split.clusterProfiler.results}
\alias{split.clusterProfiler.results}
\title{Split clusterProfiler output into one line per gene}
\usage{
\method{split}{clusterProfiler.results}(x)
}
\arguments{
\item{x}{Results from clusterProfiler. It is expected that the
clusterProfiler enrichment function was called with "readable=TRUE"}
}
\value{
data.frame with genes one per line instead of "/" separated in one
line. The rest of the original line is repeated for each gene.
}
\description{
Split clusterProfiler output into one line per gene
}
