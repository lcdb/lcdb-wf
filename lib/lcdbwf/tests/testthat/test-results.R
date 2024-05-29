library(testthat)
library(DESeq2)
library(lcdbwf)
library(rlang)
library(stringr)
library(BiocParallel)
library(ggplot2)
devtools::load_all('../../../../lib/lcdbwf')
source('test-functions.R')
config <- lcdbwf:::load_config('../../../../workflows/rnaseq/downstream/config.yaml')
text <- yaml::yaml.load_file('text.yaml')
register(MulticoreParam(config$parallel$cores))

# Mock function to capture mdcat output
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

# ------ Test build_results_tabs function ------ #
# Create objects for testing defaults
dds_and_res <- make_deseq_results()
dds <- dds_and_res$dds
res <- dds_and_res$res$res

dds_list <- list(dds1=dds)
res_list <- list(res1=list(res=res, dds='dds1', label='Defaults'))

test_that("build_results_tabs works with default config", {
  expect_silent(build_results_tabs(res_list, dds_list, config, text))
})

# Create objects for testing 'LRT'
dds_and_res <- make_deseq_results(test='LRT', reduced_design=~1)
dds <- dds_and_res$dds
res <- dds_and_res$res$res

dds_list <- list(dds1=dds)
res_list <- list(res1=list(res=res, dds='dds1', label='LRT'))

# Test build_results_tabs function
test_that("build_results_tabs works with LRT config", {
  expect_silent(build_results_tabs(res_list, dds_list, config, text))
})

test_that("build_results_tabs works with diagnostics disabled", {
  config$toggle$results_diagnostics <- FALSE
  expect_silent(build_results_tabs(res_list, dds_list, config, text))
})

test_that("build_results_tabs works with specific diagnostics results names", {
  config$toggle$results_diagnostics <- TRUE
  config$plotting$diagnostics_results_names <- c("res1")
  expect_silent(build_results_tabs(res_list, dds_list, config, text))
})

test_that("build_results_tabs handles empty res_list", {
  expect_silent(build_results_tabs(list(), dds_list, config, text))
})

test_that("check_LRT identifies LRT results correctly", {
  res_LRT_result <- create_deseq_results(test='LRT', reduced_design=~1)
  res_LRT <- res_LRT_result$res
  expect_true(check_LRT(res_LRT))

  res_Wald_result <- create_deseq_results(test='Wald')
  res_Wald <- res_Wald_result$res
  expect_false(check_LRT(res_Wald))
})

# Test that mdcat is called with expected values for LRT
test_that("build_results_tabs calls mdcat with expected values for LRT", {
  res_LRT_result <- create_deseq_results(test='LRT', reduced_design=~1)
  res_LRT <- res_LRT_result$res
  dds <- res_LRT_result$dds
  dds_list_LRT <- list(dds1=dds)
  res_list_LRT <- list(res1=list(res=res_LRT, dds='dds1', label='LRT Test Label'))

  # Capture mdcat output
  mdcat_output <<- c()
  with_mock(
    `lcdbwf:::mdcat` = mock_mdcat,
    build_results_tabs(res_list_LRT, dds_list_LRT, config, text)
  )

  expect_false(any(grepl("Wald", mdcat_output)))
  expect_true(any(grepl("LRT", mdcat_output)))
})

