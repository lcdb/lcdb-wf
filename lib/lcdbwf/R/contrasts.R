#' Helper function to get dds object from global env.
#'
#' @param dds If string, then look for that name in `dds_list` in the global
#' env. Otherwise, assume it's already a dds and immmediately return it.
get_dds <- function(dds){
  if (class(dds) != 'character'){
    return(dds)
  }

  if (!'dds_list' %in% ls(.GlobalEnv)){
    stop("Can't find dds_list in global environment.")
  }
  dds_list_local <- get("dds_list", .GlobalEnv)
  if (!(dds %in% names(dds_list_local))){
    stop(paste("Can't find name", dds, "in dds_list"))
  }
  dds <- dds_list_local[[dds]]
  return(dds)
}

#' Prepare coefficients that can be combined for more complex models.
#'
#' Heavily influenced by
#' https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md
#'
#' @param dds DESeqDataSet object or string. If string, looks in the global
#'   environment for "dds_list" and takes the dds object with this name in the
#'   list.
#' @param ... Additional filtering criteria which are handed off to
#'   dplyr::filter and used to filter the colData(dds).
#'
#' @return A named vector, representing the coefficients of the samples
#'   selected via `...`
#'
#'
#' @details
#'
#' First we simulate some data to work with:
#'
#'   dds <- makeExampleDESeqDataSet(n = 1000, m = 9, betaSD = 2)
#'   dds$condition <- NULL
#'   dds$colour <- factor(rep(c("pink", "yellow", "white"), each = 3))
#'   dds$colour <- relevel(dds$colour, "white")
#'   design(dds) <- ~colour
#'
#' Here's what the colData looks like:
#' colData(dds)
#' ## DataFrame with 9 rows and 1 column
#' ##           colour
#' ##         <factor>
#' ## sample1   pink
#' ## sample2   pink
#' ## sample3   pink
#' ## sample4   yellow
#' ## sample5   yellow
#' ## sample6   yellow
#' ## sample7   white
#' ## sample8   white
#' ## sample9   white
#'
#' If we wanted a contrast for "pigmented vs white", we could then do:
#'
#'   pigmented <- dds_coefs(dds, colour %in% c('pink', 'yellow'))
#'   white <- dds_coefs(dds, colour=='white')
#'   res <- results(dds, contrast=pigmented-white)
#'
dds_coefs <- function(dds, ..., expand=FALSE){

  dds <- get_dds(dds)

  # NOTE: throughout this function we convert to data.frame to make it easier
  # to use dplyr functions

  ci <- SummarizedExperiment::colData(dds) %>% as.data.frame()
  if (!("samplename" %in% colnames(ci))){
    stop("Need to have 'samplename' as a column in colData")
  }

  design <- DESeq2::design(dds)
  if (expand){
    mat <- DESeq2:::makeExpandedModelMatrix(dds)
  } else {
    mat <- model.matrix(design, ci)
  }

  # We will be doing some joining below but we need to know which columns to
  # select back out later
  columns <- dimnames(mat)[[2]]

  # This joins the model.matrix with the colData, applies the filter provided
  # in the ... arguments, selects back out just the model.matrix columns, and
  # returns the colMeans of the resulting filtered model.matrix.
  z <- mat %>% as.data.frame() %>%
    tibble::rownames_to_column("samplename") %>%
    dplyr::inner_join(ci, by="samplename") %>%
    dplyr::filter(...) %>%
    dplyr::select(all_of(columns)) %>%
    colMeans
  return(z)
}


#' Convenience function for building contrasts
#'
#' @description
#'
#' This function first runs DESeq2::results() and then runs DESeq2::lfcShrink()
#' on that results object. Arguments in ... are passed through to results() and
#' then to lfcShrink(), using formals() to detect which function accepts which
#' arguments. Those arguments accepted by both (contrast, lfcThreshold, format,
#' saveCols, parallel, BPPARAM) are passed to both.
#'
#' This function is useful for ensuring that the dds object used when creating
#' the results is the same one that is labeled, so that when the dds_list and
#' res_list are composed together things match up nicely.
#'
#' This function also adds the shrinkage type as an additional item to the
#' metadata of the results object so it can be reported.
#'
#' NOTE: this function expects `dds_list` and `config` objects to be available
#' in the global environment.
#'
#' @param dds_name String key into dds_list. Will be used to extract the actual
#'   dds object to provide to results()
#' @param label Label to describe this contrast which will be used in headings.
#' @param dds_list List of dds objects. If NULL, then look in the global
#'   environment for an object called "dds_list" and use that.
#' @param type Type of shrinkage for use by lfcShrink(). If no type is given,
#'   we use the current DESeq2 default argument for Type. If
#'   NULL is given, we skip lfcShrink() altogether and directly return the object from results().
#' @param ... Additional arguments are passed to results() and lfcShrink(). If
#'   "parallel" is not explicitly specified here, then look in the global env for
#'   a variable called "config" and find the parallel config setting from there.
#'
#' @return List as expected by lcdbwf, with the following names:
#'   - "res" (the results object after results() and lfcShrink())
#'   - "dds" (the string value passed in used to retrieve the actual dds object
#'     used for the results)
#'   - "label" (exactly as passed in)
make_results <- function(dds_name, label, dds_list=NULL, ...){

  if (!is.null(dds_list)){
    dds <- dds_list[[dds_name]]
  } else {
    # Get the dds object from the global dds_list.
    dds <- get_dds(dds_name)
  }

  # Save a copy of the arbitrary arguments
  dots <- list(...)

  # If not provided, detect parallel based on global config.
  if (!("parallel" %in% names(dots))){
    parallel <- FALSE
    config <- get_config()
    parallel <- config$parallel$parallel
    if (is.null(parallel)) parallel <- FALSE
  }

  # Modify the args based on what we detected.
  dots[['parallel']] <- parallel

  # Note that results() expects the argument to be called 'object' rather than,
  # say, 'dds'.
  dots[['object']] <- dds

  # Ensure any provided `test` argument is consistent with the dds object provided.
  # This uses names from mcols(dds) to detect how the dds object was created.
  test_detected <- FALSE
  if ('test' %in% names(dots)) {
    if ((dots$test == 'Wald' && any(grepl('LRT', names(S4Vectors::mcols(dds))))) ||
        (dots$test == 'LRT' && any(grepl('Wald', names(S4Vectors::mcols(dds)))))) {
      stop("The 'test' passed to make_results does not match the detected test type in dds")
    }
  } else {
    if (any(grepl('LRT', names(S4Vectors::mcols(dds))))) {
      dots$test <- 'LRT'
      test_detected <- TRUE
    } else if (any(grepl('Wald', names(S4Vectors::mcols(dds))))) {
      dots$test <- 'Wald'
      test_detected <- TRUE
    } else {
      stop("test type was missing from make_results call and could not be detected from dds")
    }
  }

  # Set the current default for 'type' from DESeq2 for lfcShrink if 'type' was not provided.
  # This inspects the function definition of lfcShrink to see what the current default is
  # (we have have seen it change before, hence the check).
  if (!'type' %in% names(dots)) {
    dots$type <- eval(formals(DESeq2::lfcShrink)$type)[1]
  }

  # Call results() with the subset of dots that it accepts.
  results_dots <- lcdbwf:::match_from_dots(dots, DESeq2::results)
  res <- do.call(DESeq2::results, results_dots)

  # When make_results is called with 'test' set to 'LRT',
  # or when make_results is called with 'test' missing but
  # DDS object contains the LRT, we convert all values in the log2FoldChange
  # column of the DESeqResults object to 0. LFC values only make sense to report for a single
  # comparison of two sample groups. This applies to the Wald test only.
  # LRT is instead performing a test of the removal of one or more factor(s) from the design formula.
  # DESeq2 reports log2FoldChange values for a single pair-wise comparison when test == 'LRT'. This
  # can be misleading and so this is our solution.

  # Adjust log2FoldChange for LRT test
  if (!is.null(dots$test) && dots$test == 'LRT') {
    res$log2FoldChange <- 0
  }

  # Checks for LRT test and non-NULL type
  if (!is.null(dots$type) && !is.null(dots$test) && dots$test == 'LRT' && !test_detected) {
    stop("You cannot pass a non-NULL or missing type to make_results with test == 'LRT'. For LRT, LFC values are set to 0 and should not be passed to lfcShrink. Use type == NULL in make_results for LRT DDS objects.")
  } else if (!is.null(dots$type) && !is.null(dots$test) && dots$test == 'LRT' && test_detected) {
    stop("You cannot pass a non-NULL or missing type to make_results with an LRT dds object. For LRT, LFC values are set to 0 and should not be passed to lfcShrink. Use type == NULL in make_results for LRT DDS objects.")
  }

  # Call lfcShrink if applicable
  if (!is.null(dots$type) && dots$test != 'LRT') {
    dots[['res']] <- res
    dots[['dds']] <- dds

    lfcShrink_dots <- lcdbwf:::match_from_dots(dots, DESeq2::lfcShrink)
    res <- do.call(DESeq2::lfcShrink, lfcShrink_dots)

    S4Vectors::metadata(res)$type <- dots$type
  } else {
    S4Vectors::metadata(res)$type <- NULL
  }

  return(
    list(
      res = res,
      dds = dds_name,
      label = label
    )
  )
}

results_diagnostics <- function(res, dds, name, config, text){
    lcdbwf:::mdcat('### Other diagnostics')
    print(knitr::kable(lcdbwf:::my_summary(res, dds, name)))

    lcdbwf:::folded_markdown(text$results_diagnostics$filter_ma, "Help")
    filterThreshold <- S4Vectors::metadata(res)$filterThreshold
    p1 <- ggplot2::ggplot(res %>% as.data.frame() %>% dplyr::mutate(filtered=res$baseMean < filterThreshold)) +
      ggplot2::aes(x=log10(baseMean), y=log2FoldChange, color=filtered) +
      ggplot2::geom_point()
    print(p1)

    lcdbwf:::folded_markdown(text$results_diagnostics$outlier_ma, "Help")
    p2 <- ggplot2::ggplot(res %>% as.data.frame() %>% dplyr::mutate(outlier=is.na(res$pvalue))) +
      ggplot2::aes(x=log10(baseMean), y=log2FoldChange, color=outlier) +
      ggplot2::geom_point()
    print(p2)

    lcdbwf:::folded_markdown(text$results_diagnostics$lfcse_basemean, "Help")
    p3 <- ggplot2::ggplot(res %>% as.data.frame() %>% dplyr::mutate(outlier=is.na(res$pvalue))) +
      ggplot2::aes(x=log10(baseMean), y=lfcSE, color=outlier) +
      ggplot2::geom_point()
    print(p3)

    lcdbwf:::folded_markdown(text$results_diagnostics$lfcse_lfc, "Help")
    p4 <- ggplot2::ggplot(res %>% as.data.frame() %>% dplyr::mutate(outlier=is.na(res$pvalue))) +
      ggplot2::aes(x=log2FoldChange, y=lfcSE, color=outlier) +
      ggplot2::geom_point()
    print(p4)

    # Save plots to a list and return for testing
    plots <- list(p1=p1, p2=p2, p3=p3, p4=p4)
    return(plots)
}
