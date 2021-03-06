% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{get.orgdb}
\alias{get.orgdb}
\title{Get the OrgDb for the specified organism, using the cached AnnotationHub.}
\usage{
get.orgdb(species, cache, annotation_key_override = NA)
}
\arguments{
\item{species}{Case-sensitive genus and species}

\item{cache}{Directory in which the AnnotationHub cache is stored}

\item{annotation_key_override}{If not NA, forces the hub to use this
accession. Use this when you know exactly which OrgDb to use.}
}
\value{
OrgDb object
}
\description{
Get the OrgDb for the specified organism, using the cached AnnotationHub.
}
