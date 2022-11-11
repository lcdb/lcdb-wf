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

  # Modify the args based on what we detected. Note that results() expects the
  dots['parallel'] <- parallel

  # Note that results() expects the argument to be called 'object' rather than,
  # say, 'dds'.
  dots[['object']] <- dds

  # Call results() with the subset of dots that it accepts.
  results_dots <- lcdbwf:::match_from_dots(dots, results)
  res <- do.call("results", results_dots)

  # We're about to call lfcShrink, but it needs the res object...so inject the
  # one we just made into dots.
  dots[['res']] <- res

  # lfcShrink also needs the dds object, so inject that too
  dots[['dds']] <- dds

  lfcShrink_dots <- lcdbwf:::match_from_dots(dots, lfcShrink)
    res <- do.call("lfcShrink", lfcShrink_dots)

  # Add the shrinkage type to the metadata of the results object.
  #
  # If "type" was specified when calling this function, it's easy and we use
  # that. Otherwise, if it was not specified then DESeq2 used the default.
  # Since that default can change as we have seen in the past, we need to
  # inspect the lfcShrink function itself to see what the current default is,
  # and use that.
  shrinkage_type <- dots[['type']]
  if (is.null(shrinkage_type)){
    # The definition of lfcShrink has a character vector as the type argument,
    # and we want to extract the first thing in that vector. But formals()
    # return strings, so we need to eval that string to convert it to
    # a character vector such that we can extract the first thing.
    #
    # In recent versions this should evaluate to "apeglm". But this way we
    # are inspecting the function itself if it ever changes.
    shrinkage_type <- eval(formals(DESeq2::lfcShrink)$type)[1]
  }

  # Add to results object so we can report it out later.
  metadata(res)$type <- shrinkage_type

  return(
    list(
      res=res,
      dds=dds_name,
      label=label
    )
  )
}


results_diagnostics <- function(res, dds, name, config, text){
    lcdbwf:::mdcat('### Other diagnostics')
    print(knitr::kable(lcdbwf:::my_summary(res, dds, name)))

    lcdbwf:::folded_markdown(text$results_diagnostics$filter_ma, "Help")
    filterThreshold <- metadata(res)$filterThreshold
    p <- ggplot(res %>% as.data.frame() %>% mutate(filtered=res$baseMean < filterThreshold)) +
      aes(x=log10(baseMean), y=log2FoldChange, color=filtered) +
      geom_point()
    print(p)

    lcdbwf:::folded_markdown(text$results_diagnostics$outlier_ma, "Help")
    p <- ggplot(res %>% as.data.frame() %>% mutate(outlier=is.na(res$pvalue))) +
      aes(x=log10(baseMean), y=log2FoldChange, color=outlier) +
      geom_point()
    print(p)

    lcdbwf:::folded_markdown(text$results_diagnostics$lfcse_basemean, "Help")
    p <- ggplot(res %>% as.data.frame() %>% mutate(outlier=is.na(res$pvalue))) +
      aes(x=log10(baseMean), y=lfcSE, color=outlier) +
      geom_point()
    print(p)

    lcdbwf:::folded_markdown(text$results_diagnostics$lfcse_lfc, "Help")
    p <- ggplot(res %>% as.data.frame() %>% mutate(outlier=is.na(res$pvalue))) +
      aes(x=log2FoldChange, y=lfcSE, color=outlier) +
      geom_point()
    print(p)
}
