#' Build tabbed output for results
#'
#' This function outputs the RMarkdown and plots for building nested tabs for
#' each set of results. Make sure this is called within a `results="asis"`
#' chunk.
#'
#' @param res_list List of results objects
#' @param dds_list List of dds objects
#' @param config Config object. Use config$plotting$diagnostics_results_names
#'   to configure which results names (default is all). Use
#'   config$toggle$results_diagnostics to enable/disable diagnostics.
#' @param text Text config object. text$results_plots is the relevant section used by this function.
#'
#' @return NULL, but creates the plots and tabs as a side effect.
build_results_tabs <- function(res_list, dds_list, config, text){

  # Default to no diagnostics
  diagnostics_names <- c()

  # If diagnostics are toggled on, then default running them on all results
  # names. Otherwise restrict to what was provided.
  if (config$toggle$results_diagnostics){
    if (
        is.null(config$plotting$diagnostics_results_names) ||
          length(config$plotting$diagnostics_results_names) == 0){
      diagnostics_names <- names(res_list)
    } else {
      diagnostics_names <- config$plotting$diagnostics_results_names
    }
  }

  for (name in names(res_list)){
    dds_i <- dds_list[[res_list[[name]][['dds']] ]]
    res_i <- res_list[[name]][['res']]
    label <- res_list[[name]][['label']]
    # Do any contrasts contain LRT?
    contains_LRT <- check_LRT(res_i)
    genes_to_label <- lcdbwf:::genes_to_label(res_i, n=5, config)
    lcdbwf:::mdcat('## ', label, ' {.tabset}')

    # TODO:
    # https://hdgitappip01.nichd.nih.gov/bspc/core/labeled-ma-plot/-/blob/master/plotMA.label.R
    # is the latest version

    lcdbwf:::mdcat('### M-A plot')
    lcdbwf:::folded_markdown(text$results_plots$ma, "Help")
    # If any contrasts contain LRT, print the source of LFC
    # and p values above MA & Volcano plots
    if (contains_LRT) {
      mdcat("LRT log2FoldChange values have been set to 0")
    }
    print(lcdbwf:::plotMA_label(
      res_i,
      genes_to_label=genes_to_label,
      label_column=config$annotation$label_column))

    lcdbwf:::mdcat('### Volcano plot')
    lcdbwf:::folded_markdown(text$results_plots$volcano, "Help")
    if (contains_LRT) {
      mdcat("LRT log2FoldChange values have been set to 0")
    }

    print(lcdbwf:::plot_volcano_label(
      res_i,
      genes_to_label=genes_to_label,
      label_column=config$annotation$label_column))

    lcdbwf:::mdcat('### P-value distribution')
    lcdbwf:::folded_markdown(text$results_plots$pval_hist, "Help")
    print(lcdbwf:::pval_hist(res_i))

    if (config$toggle$results_diagnostics){
      lcdbwf:::results_diagnostics(res=res_i, dds=res_list[[name]]$dds, name=name, config=config, text=text)
    }
  }
}

#' Check for LRT in a results object's metadata
#' @param res_i DESeqResults object
#' @return Boolean TRUE if results object's pvalues were provided by
#' nBinomLRT and FALSE if the Wald test was used.
check_LRT <- function(res_i) {
  mcols_pval <- mcols(res_i)$description[9]
  mcols_pval <- grepl('LRT', mcols_pval)
  return(mcols_pval)
}
