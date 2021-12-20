#' Build tabbed output for results
#'
#' This function outputs the RMarkdown and plots for building nested tabs for
#' each set of results. Make sure this is called within a `results="asis"`
#' chunk.
#'
#' @param res_list List of results objects
#' @param dds_list List of dds objects
#' @param config Config object
#' @param text Text config object
#'
#' @return NULL, but creates the plots and tabs as a side effect.
build_results_tabs <- function(res_list, dds_list, config, text){

  # Default to no diagnostics
  diagnostics_names <- c()

  # If not oth
  if (config$toggle$results_diagnostics){
    if (
        is.null(config$plotting$diagnostics_results_names) ||
          length(config$plotting$diagnostics_results_names) == 0){
      diagnostics_names <- names(res_list)
    } else {
      diagnostics_names <- config$plotting$diagnostics_results_names
  }

  for (name in names(res_list)){
    dds_i <- dds_list[[res_list[[name]][['dds']] ]]
    res_i <- res_list[[name]][['res']]
    label <- res_list[[name]][['label']]
    genes_to_label <- lcdbwf::genes_to_label(res_i, n=5, config)
    lcdbwf::mdcat('## ', label, ' {.tabset}')

    # TODO:
    # https://hdgitappip01.nichd.nih.gov/bspc/core/labeled-ma-plot/-/blob/master/plotMA.label.R
    # is the latest version

    lcdbwf::mdcat('### M-A plot')
    lcdbwf::folded_markdown(text$results_plots$ma, "Help")
    print(lcdbwf::plotMA_label(
      res_i,
      genes_to_label=genes_to_label,
      label_column=config$annotation$label_column))

    lcdbwf::mdcat('### Volcano plot')
    lcdbwf::folded_markdown(text$results_plots$volcano, "Help")
    print(lcdbwf::plot_volcano_label(
      res_i,
      genes_to_label=genes_to_label,
      label_column=config$annotation$label_column))

    lcdbwf::mdcat('### P-value distribution')
    lcdbwf::folded_markdown(text$results_plots$pval_hist, "Help")
    lcdbwf::pval_hist(res_i)

    if (config$toggle$results_diagnostics){
      lcdbwf::results_diagnostics(res=res_i, dds=res_list[[name]]$dds, name=name, config=config, text=text)
    }
  }
}
