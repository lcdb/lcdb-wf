#' Helper function to get dds object from global env.
#'
#' @param dds If string, then look for that name in `dds_list` in the global
#'   env. Otherwise, return it as-is.
get_dds <- function(dds){
  if (class(dds) == 'character'){
    if (!'dds_list' %in% ls(.GlobalEnv)){
      stop("Can't find dds_list in global environment.")
    }
    dds_list_local <- get("dds_list", .GlobalEnv)

    if (!(dds %in% names(dds_list_local))){
      stop(paste("Can't find name", dds, "in dds_list"))
    }
    dds <- dds_list_local[[dds]]
  }
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


#' Convenience function for building contrasts using the data structures
#' expected by lcdbwf.
#'
#' NOTE: this function makes some assumptions about the global environment. For
#' example, it looks for the config object there, and expects a list called
#' "dds_list" to exist as well. Currently the config is only used for
#' parallelization.
#'
#' This function is useful for ensuring that the dds object used when creating
#' the results is the same one that is labeled, so that when the dds_list and
#' res_list are composed together things match up nicely.
#'
#' @param dds_name String key into dds_list
#' @param contrast Any contrast supported by DESeq2::lfcShrink
#' @param label Label or description to describe this contrast
#' @param ... Additional arguments are passed to lfcShrink
#'
#' @return List with names "res", "dds", and "label".
make_results <- function(dds_name, label, ...){

  dds <- get_dds(dds_name)


  # Modify the args based on what we detected
  dots <- list(...)
  dots[['object']] <- dds

  if (!"parallel" %in% names(dots)){
    parallel <- FALSE
    config <- get_config()
    parallel <- config$parallel$parallel
    if (is.null(parallel)) parallel <- FALSE
  }

  #if (parallel) dots[['parallel']] <- parallel

  # Call results()
  res_dots <- lcdbwf::match_from_dots(dots, results)
  res <- do.call("results", lcdbwf::match_from_dots(dots, results))

  # We're about to call lfcShrink, but it needs the res object...so inject the
  # one we just made into dots.
  dots[['res']] <- res
  dots[['dds']] <- dds
  res <- do.call("lfcShrink", lcdbwf::match_from_dots(dots, lfcShrink))


  return(
    list(
      res=res,
      dds=dds_name,
      label=label
    )
  )
}

