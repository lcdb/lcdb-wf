# ------------------------------------------------------------------------------
# Functions for working with DESeqDataSet (dds) objects
# ------------------------------------------------------------------------------


salmon.path.func <- function (x) file.path('..', 'data', 'rnaseq_samples', x, paste0(x, '.salmon'), 'quant.sf')
kallisto.path.func <- function (x) file.path('..', 'data', 'rnaseq_samples', x, paste0(x, '.salmon'), 'quant.sf')



#' Make a single dds object
#'
#' This function is intended to be applied to a list of design data along with
#' a config object to create a list of dds objects:
#'
#'    dds_list <- map(dds_list, make_dds, config)
#'
#' While the intended use is to provide a config object, individual arguments
#' can be used to override individual config items. For these arguments, the
#' default in the function signature is left as NULL (even if the value should
#' be TRUE or FALSE) so as to detect whether or not the config should be
#' overridden.
#'
#' @param design_data Named list of 2 to 4 items. The names "sampletable" and
#' "design" are required.  The optional named items are "file", which is the
#' featureCounts file containing counts for all samples, and "args" which is
#' a list of arguments to be passed to the constructor (e.g.,
#' `args=list(subset.counts=TRUE))`.
#'
#' @param collapse_by The column to collapse technical replicates by. Rows in
#'   the sampletable that share the same value of this column will be combined
#'   using DESeq2::collapseReplicates.
#'
#' @param strip_dotted_version If TRUE, then remove Ensembl-style dotted
#'   version numbers from gene IDs (ENSG000001.1 -> ENSG000001)
#'
#' @param salmon_pattern, kallisto_pattern Specify the patterns to locations of
#'   Salmon or Kallisto files. Use the special placeholder string
#'   `__SAMPLENAME__` which will be replaced with the sample name. Only
#'   relevant if one of config$toggle$salmon or config$toggle$kallisto are
#'   TRUE.
#'
#' @param ... Additional arguments will be passed on to the DESeq() call (e.g.,
#'   parallel, fitType, etc)
#'
#' @param featureCounts Location of featureCounts output to be loaded
make_dds <- function(design_data, config=NULL, collapse_by=NULL,
                     strip_dotted_version=NULL,
                     featureCounts='../data/rnaseq_aggregation/featurecounts.txt',
                     salmon_pattern="../data/rnaseq_samples/__SAMPLENAME__/__SAMPLENAME__.salmon/quant.sf",
                     kallisto_pattern="../data/rnaseq_samples/__SAMPLENAME__/__SAMPLENAME__.kallisto/abundance.h5",
                     ...){

  # Note we're using pluck() here for the conveneience of setting defaults
  coldata <- pluck(design_data, 'sampletable')
  design <- pluck(design_data, 'design')
  location <- pluck(design_data, 'file', .default=featureCounts)
  salmon <- pluck(design_data, 'salmon')
  kallisto <- pluck(design_data, 'kallisto')
  arg_list <- pluck(design_data, 'args')

  # Allow overriding of config values.
  if (!is.null(config)){
    if (is.null(salmon)) salmon <- config$toggle$salmon
    if (is.null(kallisto)) kallisto <- config$toggle$kallisto
    if (is.null(collapse_by)) collapse_by <- config$main$collapse_by
    if (is.null(strip_dotted_version)) strip_dotted_version <- config$main$strip_dotted_version
  }

  if (salmon & kallisto){
    stop("Both salmon and kallisto are set to TRUE, not sure how to handle this.")
  }

  if (salmon | kallisto){
    # If these arguments were provided, the corresponding loading functions
    # don't accept them so we need to remove. Issue a warning as well.
    if (!is.null(arg_list$subset.counts) | !is.null(arg_list$sample.func)){
      warning("Salmon or Kallisto was specified, but additional arguments ",
              "were provided to the loading function.")
      arg_list$subset.counts <- NULL
      arg_list$sample.func <- NULL
    }

    # Next, we need a tx2gene dataframe. We can get this from a TxDb, which in
    # turn can be retrieved from AnnotationHub, which in turn can be configured
    # with the config object. Luckily, we have it here!
    txdb <- get_annotation_db(config, dbtype="TxDb")
    k <- keys(txdb, keytype="TXNAME")
    tx2gene <- select(txdb, k, "GENEID", "TXNAME")

  }

  if (salmon){
      coldata$salmon.path <- sapply(coldata$samplename, function (x) gsub("__SAMPLENAME__", x, salmon_pattern))
      txi <- tximport(coldata[, 'salmon.path'], type='salmon', tx2gene=tx2gene, ignoreTxVersion=strip_dotted_version)
      dds <- DESeqDataSetFromTximport(txi, colData=coldata, design=design)

  } else if (kallisto) {
      coldata$kallisto.path <- sapply(coldata$samplename, function (x) gsub("__SAMPLENAME__", x, kallisto_pattern))
      txi <- tximport(coldata[, 'kallisto.path'], type='kallisto', tx2gene=tx2gene, ignoreTxVersion=strip_dotted_version)
      dds <- DESeqDataSetFromTximport(txi, colData=coldata, design=design)

  } else {
      dds <- exec(
          DESeqDataSetFromCombinedFeatureCounts,
              location,
              sampletable=coldata,
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
