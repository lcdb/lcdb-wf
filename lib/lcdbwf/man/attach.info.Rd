% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{attach.info}
\alias{attach.info}
\title{Attach additional information from an OrgDb to a results dataframe}
\usage{
attach.info(
  res,
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "UNIPROT", "ALIAS")
)
}
\arguments{
\item{res}{DESeq2 results object}

\item{keytype}{Key type for rownames(res)}

\item{columns}{Character vector of columns to add from OrgDb to results}
}
\description{
Attach additional information from an OrgDb to a results dataframe
}
