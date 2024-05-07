#' Load a single combined featureCounts table into a DESeq object.
#'
#' @param filename Filename containing featureCounts combined output
#' @param sampletable data.frame containing at least sample names as first
#'        column
#' @param design Model used for creating DESeq object
#' @param sample_func Function that will be applied to each column name
#'        to align it to the input sampletable. Takes a character vector as
#'        input and returns a character vector of adjusted strings.
#' @param subset_counts If TRUE, then counts are subsetted by the provided
#'        sampletable. That is, only the counts columns that match a rowname of
#'        the provided sampletable are included. If FALSE (default), an error
#'        is raised reporting the differences in sampletable and counts data.
#'        It is always an error if counts data does not have an entry for
#'        a sample in the sampletable.
#'
#' @return DESeq object
#'
#' Additional args are passed to DESeq2::DESeqDataSetFromMatrix.
DESeqDataSetFromCombinedFeatureCounts <- function(filename, sampletable,
                                                design,
                                                sample_func=lcdbwf_samplename,
                                                subset_counts=FALSE, ...){

  if (is.null(subset_counts)){
    subset_counts <- FALSE
  }
  # The sampletable may be data.frame or tibble; if it's a tibble then it
  # likely doesn't have rownames. So in this function we assume that it's the
  # first column that contains the samplenames.

  # Read in the counts TSV, use gene ID as rownames, and get rid of the Chr,
  # Start, End, Strand, Length columns
  m <- readr::read_tsv(filename, comment="#", col_types=readr::cols()) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames('Geneid') %>%
      dplyr::select(-(1:5)) %>%
      as.data.frame()

  # The column names of the imported table are not guaranteed to match the
  # samples in the metadata, so we need to make sure things match up
  # correctly when renaming columns.
  #
  # sample_func should convert filenames (which are columns in featurecounts
  # table) with samples (as listed in the sampletable)
  x <- colnames(m) %>% sample_func()

  samplenames <- sampletable[,1]
  counts.not.sampletable <- setdiff(x, samplenames)
  sampletable.not.counts <- setdiff(samplenames, x)


  if (!all(x %in% samplenames)){
    # sampletable is missing:
    if (!subset_counts){
      stop(
        paste(
          'The following samples are in the counts data but not the',
          'sampletable. If this is intended, consider using `subset_counts=TRUE`',
          'to remove them from the counts:', paste(counts.not.sampletable, collapse=', ')
        )
      )
    }
  }
  if (!all(samplenames %in% x)){
    stop(
     paste(
       'The following samples are in the sampletable but not in the counts data. Check sample_func?',
       paste(sampletable.not.counts, collapse=', '))
    )
  }

  colnames(m) <- x

  # We should be good -- so reorder columns in `m` to match sampletable
  m <- m[, sampletable[,1]]

  # Before creating the dds object, we will drop levels from factors, in case
  # that was not done before sending in the (possibly subsetted) sampletable.

  sampletable <- droplevels(as.data.frame(sampletable))

  object <- DESeq2::DESeqDataSetFromMatrix(countData=m, colData=sampletable, design=design, ...)
  return(object)
}

#' Load featureCounts output into a DESeq object.
#'
#' Revised version of DESeq2::DESeqDataSetFromHTSeqCount to handle the
#' featureCounts default output format, which contains many more columns.
#'
#' @param sampleTable data.frame containing at least "featurecounts.path" column
#' @param directory Paths to featureCounts output are relative to this dir
#' @param design Model used for creating DESeq object
#'
#' @return DESeq object
#'
#' Additional args are passed to DESeq2::DESeqDataSetFromMatrix.
DESeqDataSetFromFeatureCounts <- function (sampleTable, directory='.', design,
                                           ignoreRank=FALSE,  ...)
{
  l <- lapply(
    as.character(sampleTable[, 'featurecounts.path']),
    function(fn) read.table(file.path(directory, fn), stringsAsFactors=FALSE, skip=2)
  )
  if (!all(sapply(l, function(a) all(a$V1 == l[[1]]$V1))))
    stop("Gene IDs in first column differ between files")
  tbl <- sapply(l, function(a) a$V7)
  colnames(tbl) <- sampleTable[, 1]
  rownames(tbl) <- l[[1]]$V1
  rownames(sampleTable) <- sampleTable[, 1]
  object <- DESeq2::DESeqDataSetFromMatrix(
    countData=tbl,
    colData=sampleTable[, -grepl('path', colnames(sampleTable)), drop=FALSE],
    design=design,
    ignoreRank,
    ...)
  return(object)
}
