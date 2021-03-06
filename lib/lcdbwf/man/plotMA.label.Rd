% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{plotMA.label}
\alias{plotMA.label}
\title{Plot a MA plot labeled with selected genes}
\usage{
plotMA.label(
  res,
  fdr.thres = 0.1,
  fc.thres = 0,
  fc.lim = NULL,
  genes.to.label = NULL,
  col = NULL
)
}
\arguments{
\item{fdr.thres}{FDR threshold for defining statistical significance}

\item{fc.thres}{log2FoldChange cutoff for defining statistical significance}

\item{fc.lim}{User-defined limits for y-axis (log2FC). If NULL, this is defined as
the (floor, ceiling) of range(res$log2FoldChange)}

\item{genes.to.label}{Genes to label on the MA plot. NULL, by default}

\item{col}{Column of \code{res} in which to look for \code{genes.to.label}. If NULL,
rownames are used.}

\item{res.list}{data.frame (or list) with log2FoldChange and padj columns (or elements).
Row names of the data.frame should be gene IDs/symbols and match the format of lab.genes}
}
\value{
Handle to ggplot
}
\description{
Plot a MA plot labeled with selected genes
}
