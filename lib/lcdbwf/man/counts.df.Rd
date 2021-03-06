% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{counts.df}
\alias{counts.df}
\title{Make dataframe of normalized gene counts}
\usage{
counts.df(
  dds,
  res,
  sel.genes = NULL,
  label = NULL,
  rank.col = "padj",
  pc = 0.5
)
}
\arguments{
\item{dds}{DESeq2 object}

\item{res}{DESeq2 results object}

\item{sel.genes}{list of genes to consider}

\item{label}{column(s) to be included to the plot labels}

\item{pc}{count number to be added to the normalized counts
Typically a pc of 0.5 is added to allow plotting in log scale}
}
\value{
dataframe
}
\description{
Make dataframe of normalized gene counts
}
