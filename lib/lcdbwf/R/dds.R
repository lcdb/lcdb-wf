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
#' "design" are required.  The optional named items are "filename", which is the
#' featureCounts file containing counts for all samples, and "args" which is
#' a list of arguments to be passed to the constructor (e.g.,
#' `subset_counts=TRUE`).
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

  # Note we're using pluck() here for the convenience of setting defaults
  

  coldata <- purrr::pluck(design_data, 'sampletable')
  design <- purrr::pluck(design_data, 'design')
  location <- purrr::pluck(design_data, 'filename', .default=featureCounts)
  salmon <- purrr::pluck(design_data, 'salmon')
  kallisto <- purrr::pluck(design_data, 'kallisto')
  subset_counts <- purrr::pluck(design_data, 'subset_counts')
  sample_func <- purrr::pluck(design_data, 'sample_func', .default=lcdbwf_samplename)

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
    if (!is.null(subset_counts) | !is.null(sample_func)){
      warning("Salmon or Kallisto was specified, but additional arguments ",
              "were provided to the loading function.")
      subset_counts <- NULL
      sample_func <- NULL
    }

    # For Salmon and Kallisto, we need a tx2gene dataframe. We can get this
    # from a TxDb, which in turn can be retrieved from AnnotationHub, which in
    # turn can be configured with the config object. Luckily, we have the
    # config object here!
    txdb <- get_annotation_db(config, dbtype="TxDb")
    k <- keys(txdb, keytype="TXNAME")
    tx2gene <- select(txdb, k, "GENEID", "TXNAME")
  }

  if (salmon){
      coldata$salmon.path <- sapply(coldata$samplename, function (x) gsub("__SAMPLENAME__", x, salmon_pattern))
      txi <- tximport::tximport(coldata[, 'salmon.path'], type='salmon', tx2gene=tx2gene, ignoreTxVersion=strip_dotted_version)
      dds <- DESeq2::DESeqDataSetFromTximport(txi, colData=coldata, design=design)

  } else if (kallisto) {
      coldata$kallisto.path <- sapply(coldata$samplename, function (x) gsub("__SAMPLENAME__", x, kallisto_pattern))
      txi <- tximport::tximport(coldata[, 'kallisto.path'], type='kallisto', tx2gene=tx2gene, ignoreTxVersion=strip_dotted_version)
      dds <- DESeq2::DESeqDataSetFromTximport(txi, colData=coldata, design=design)

  } else {
    dds <- lcdbwf:::DESeqDataSetFromCombinedFeatureCounts(
      location,
      sampletable=coldata,
      design=design,
      sample_func=sample_func,
      subset_counts=subset_counts
    )
  }

  if (strip_dotted_version){
    dds <- lcdbwf:::strip_dotted_version_from_dds(dds)
  }

  if(!is.null(collapse_by)){
      dds <- lcdbwf:::collapseReplicates2(dds, dds[[collapse_by]])
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
strip_dotted_version_from_dds <- function(dds, force=FALSE){

  check <- head(rownames(dds), 500)
  if (!all(grepl("^ENS", check)) && !force){
    stop(paste("Gene names don't appear to be Ensembl, if you really want",
               "to remove versions, then use force=TRUE"))
  }
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


#' Plot various diagnostics for dds objects.
#'
#' @param dds_list List of dds objects
#' @param text Text config object
dds_diagnostics <- function(dds_list, text){
  mdcat('## Other diagnostics {.tabset}')
  mdcat('This section provides details on the DESeqDataSet object created above. These ',
        'can be useful for diagnosing any issues.')

  for (name in names(dds_list)){
    mdcat("### ", name, '{.tabset}')

    mdcat("#### Dispersion estimates")
    mdcat(text$dds_diagnostics$dispersion)
    DESeq2::plotDispEsts(dds_list[[name]])

    mdcat("#### Sparsity plot")
    mdcat(text$dds_diagnostics$sparsity)
    print(lcdbwf:::plotSparsity2(dds_list[[name]]))

    mdcat("#### Outliers")
    mdcat(text$dds_diagnostics$outliers)
    p <- assays(dds_list[[name]])[['cooks']] %>%
      as.data.frame() %>%
      tidyr::pivot_longer(everything()) %>%
      ggplot() +
        aes(x=name, y=log10(value)) +
        geom_boxplot() +
        ylab("log10(Cook's distance)")
    print(p)

    mdcat("#### colData")
    mdcat(text$dds_diagnostics$colData)
    cdata <- colData(dds_list[[name]]) %>% as.data.frame
    cdata <- cdata[, !grepl("filename", colnames(cdata))]
    print(htmltools::tagList(datatable(cdata)))

    mdcat("#### Design matrix")
    mdcat(text$dds_diagnostics$design_matrix, " The design is: `", deparse(design(dds_list[[name]])), "`")
    mmat <- model.matrix(design(dds_list[[name]]), data=colData(dds_list[[name]])) %>% as.data.frame()
    print(htmltools::tagList(datatable(mmat)))
  }
}
