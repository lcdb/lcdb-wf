# ------------------------------------------------------------------------------
# Functions for working with DESeqDataSet (dds) objects
# ------------------------------------------------------------------------------
#
#' Make a single dds object
#'
#' This function is intended to be applied to a list of design data along with
#' a config object to create a list of dds objects:
#'
#'    dds_list <- map(dds_list, make_dds, config)
#'
#' While the intended use is to provide a config object, individual arguments
#' can be used to override individual config items.
#'
#' @param design_data Named list of 2 to 4 items. At least "sampletable" (the
#'   colData, as a data.frame or tibble) and "design" (e.g., `~group`) are
#'   required. The optional named items are "file", which is the featureCounts
#'   file containing counts for all samples, and "args" which is a list of
#'   arguments to be passed to the constructor (e.g.,
#'   `args=list(subset.counts=TRUE))`.
#'
#' @param salmon Set to TRUE if you want to create dds objects from Salmon
#' counts (will be loaded via tximport) Otherwise, use featureCounts
#'
#' @param collapse_by The column to collapse technical replicates by. Rows in
#'   the sampletable that share the same value of this column will be combined
#'   using DESeq2::collapseReplicates.
#'
#' @param strip_dotted_version If TRUE, then remove Ensembl-style dotted
#'   version numbers from gene IDs (ENSG000001.1 -> ENSG000001)
#'
#' @param ... Additional arguments will be passed on to the DESeq() call (e.g.,
#'   parallel, fitType, etc)
#'
#' @param featureCounts Location of featureCounts output to be loaded
make_dds <- function(design_data, config=NULL, salmon=NULL, collapse_by=NULL,
                     strip_dotted_version=NULL,
                     featureCounts='../data/rnaseq_aggregation/featurecounts.txt',
                     ...){

  # Note we're using pluck() here for the conveneience of setting defaults
  colData <- pluck(design_data, 'sampletable')
  design <- pluck(design_data, 'design')
  location <- pluck(design_data, 'file', .default=featureCounts)
  arg_list <- pluck(design_data, 'args')

  # Allow overriding of config values.
  if (!is.null(config)){
    if (is.null(salmon)) salmon <- config$toggle$salmon
    if (is.null(collapse_by)) collapse_by <- config$main$collapse_by
    if (is.null(strip_dotted_version)) strip_dotted_version <- config$main$strip_dotted_version
  }

  if (salmon){
      dds <- exec(
          DESeqDataSetFromTximport,
          salmon.files,
          colData=colData[, -grepl('path', colnames(colData)), drop=FALSE],
          design=design,
          !!!arg_list)
  } else {
      dds <- exec(
          DESeqDataSetFromCombinedFeatureCounts,
              location,
              sampletable=colData,
              design=design,
              !!!arg_list)
  }

  if (strip_dotted_version){
    dds <- lcdbwf::strip_dotted_version_from_dds(dds)
  }

  if(!is.null(collapse_by)){
      dds <- lcdbwf::collapseReplicates2(dds, dds[[collapse_by]])
  }

  dds <- DESeq(dds, ...)
  return(dds)
}

#' Strip dotted version off of the rownames of a dds object
#'
#' Ensembl annotations frequently come with a dotted version number (e.g.,
#' ENSG00000001.3), but OrgDbs do not. This function removes the dotted
#' versions from the gene IDs of a DESeqDataSet object.
#'
#' @param dds DESeqDataSet object
strip_dotted_version_from_dds <- function(dds){
  rownames(dds) <- sapply(
    strsplit(rownames(dds),'.', fixed=TRUE),
    function (x) {ifelse(grepl('_', x[2]), paste(x[1], x[2], sep='.'), x[1])}
  )
    return(dds)
}


#' DESeq2::collapseReplicates, but also fix the first column
#'
#' DESeq2::collapseReplicates returns an object whose colData contains the
#' column used for collapsing, but only the first unique value of a collapsed
#' group is returned. This function makes the first column the same as the
#' rownames. This in turn allows the colData to meet expectations of other
#' lcdbwf functions and play nicer with dplyr.
#'
#' @param object Object to collapse. Typically a DESeq2 dds object
#' @param groupby Factor to group by. Typically a column from dds indicating
#'   biological replicate (e.g., dds$biorep)
collapseReplicates2 <- function(object, groupby){
    collapsed <- DESeq2::collapseReplicates(object, groupby)
    colData(collapsed)[,1] <- rownames(colData(collapsed))
    return(collapsed)
}

