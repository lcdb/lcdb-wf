source('test-functions.R')
config <- lcdbwf:::load_config('../../../../workflows/rnaseq/downstream/config.yaml')
text <- yaml::yaml.load_file('../../../../workflows/rnaseq/downstream/text.yaml')

# Wrapper function for inherits ggplot
is_ggplot <- function(x) {
  inherits(x, "ggplot")
}

# Function to capture mdcat output
mdcat_output <- c()
mock_mdcat <- function(...) {
  mdcat_output <<- c(mdcat_output, paste(..., collapse = " "))
}

# Helper function to create DESeqResults object
make_deseq_results <- function(test='Wald', type='ashr', reduced_design=NULL) {
  design_data <- make_design_data()
  design_data$test <- test
  design_data$reduced_design <- reduced_design
  label <- paste0("test=", test %||% "NULL/default (Wald)", ", type=", type %||% "NULL (Skip)")
  dds <- make_dds(design_data,
                  config=config,
                  featureCounts='featurecounts.txt',
                  parallel=config$parallel$parallel)
  tmp_dds_list = list(dds=dds)
  res <- make_results(dds_name='dds',
                      label=label,
                      dds_list=tmp_dds_list,
                      type=type)
  return(list(dds=dds, res=res))
}

# Create objects for testing defaults
dds_and_res <- make_deseq_results()
dds <- dds_and_res$dds
res <- dds_and_res$res$res

wald_dds_list <- list(dds1=dds)
wald_res_list <- list(res1=list(res=res, dds='dds1', label='Defaults'))
wald_res_list <- lcdbwf:::attach_extra(wald_res_list, config)

# Create objects for testing 'LRT'
dds_and_res <- make_deseq_results(test='LRT', type=NULL, reduced_design=~1)
dds <- dds_and_res$dds
res <- dds_and_res$res$res

lrt_dds_list <- list(dds1=dds)
lrt_res_list <- list(res1=list(res=res, dds='dds1', label='LRT'))
lrt_res_list <- lcdbwf:::attach_extra(lrt_res_list, config)

# ------ Test build_results_tabs function ------ #
test_that("build_results_tabs works with Wald test", {
  # build_results_tabs requires 'dds_list' in .GlobalEnv
  dds_list <<- wald_dds_list
  plots <- build_results_tabs(wald_res_list, wald_dds_list, config, text)

  # Check that each plot in the list is a ggplot object
  for (name in names(plots)) {
    expect_true(is_ggplot(plots[[name]]$ma_plot))
    expect_true(is_ggplot(plots[[name]]$volcano_plot))
    expect_true(is_ggplot(plots[[name]]$pval_hist_plot))
    # Check diagnostic plots
    if (config$toggle$results_diagnostics) {
      # diag_plot_list is a list of ggplot objects
      for (diag_plot in plots[[name]]$diag_plot_list) {
        expect_true(is_ggplot(diag_plot))
      } # for diag_plot
    } # if config
  } # for name
}) # test_that

test_that("build_results_tabs works with LRT", {
  # build_results_tabs requires 'dds_list' in .GlobalEnv
  dds_list <<- lrt_dds_list
  plots <- build_results_tabs(lrt_res_list, lrt_dds_list, config, text)

  # Check that each plot in the list is a ggplot object
  for (name in names(plots)) {
    expect_true(is_ggplot(plots[[name]]$ma_plot))
    expect_true(is_ggplot(plots[[name]]$volcano_plot))
    expect_true(is_ggplot(plots[[name]]$pval_hist_plot))
    # Check diagnostic plots
    if (config$toggle$results_diagnostics) {
      # diag_plot_list is a list of ggplot objects
      for (diag_plot in plots[[name]]$diag_plot_list) {
        expect_true(is_ggplot(diag_plot))
      } # for diag_plot
    } # if config
  } # for name
}) # test_that
# ---------------------------------------------- #

test_that("check_LRT identifies LRT results correctly", {
  expect_true(check_LRT(lrt_res_list$res1$res))
  expect_false(check_LRT(wald_res_list$res1$res))
})

# Test that mdcat is called with expected values for LRT
test_that("build_results_tabs calls mdcat with expected character for LRT", {
  # build_results_tabs requires 'dds_list' in .GlobalEnv
  dds_list <<- lrt_dds_list
  output <- capture_output({
    suppressWarnings(build_results_tabs(lrt_res_list, lrt_dds_list, config, text))
  })
  expect_true(any(grepl("LRT log2FoldChange values have been set to 0", output)))
})

# Test that mdcat is not called with LRT expected values for Wald test
test_that("build_results_tabs does not call mdcat with LRT expected character for Wald", {
  # build_results_tabs requires 'dds_list' in .GlobalEnv
  dds_list <<- wald_dds_list
  output <- capture_output({
    build_results_tabs(wald_res_list, wald_dds_list, config, text)
  })
  expect_false(any(grepl("LRT log2FoldChange values have been set to 0", output)))
})
