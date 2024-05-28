library(DESeq2)
library(testthat)
library(rlang)
library(stringr)
library(BiocParallel)
devtools::load_all('../../../../lib/lcdbwf')
config <- lcdbwf:::load_config('config.yaml')
source('test-functions.R')
register(MulticoreParam(config$parallel$cores))

# Test all combinations of test and type
# NULL shrinkage type skips lfcShrink
# NULL test type runs Wald (default test)
tests <- list('Wald', 'LRT', NULL)
shrinkage_types <- list('ashr', 'apeglm', 'normal', NULL)
contrast <- c("condition", "treatment", "control")
coef <- "condition_treatment_vs_control"
dds_list <- make_dds_list()
lrt_design_data <- make_lrt_design_data()

### TESTING ###
#test <- 'Wald'
#type <- 'ashr'
#dds_name <- 'dds_wald'
#contrast <- c("condition", "treatment", "control")
#label <- paste0("test=", test %||% "NULL/default (Wald)", ", type=", type %||% "NULL (Skip)")
##############

for (test in tests) {
  for (type in shrinkage_types) {
    if (test == 'Wald' || is.null(test)) {
      dds_name <- 'dds_wald'
    } else if (test == 'LRT') {
      dds_name <- 'dds_lrt'
    }
    label <- paste0("test=", test %||% "NULL/default (Wald)", ", type=", type %||% "NULL (Skip)")
    test_that(paste("make_results works correctly with", label), {
      if ((!is.null(test) && test == 'LRT') && is.null(type)) {
        res <- make_results(dds_name=dds_name, label=label, type=NULL) # No contrast when running test == 'LRT'
        check_results(res, lrt_design_data, label, test=test, type=NULL)
      } else if ((!is.null(test) && test == 'LRT') && (!is.null(type) && !type %in% c('apeglm','normal'))) {
        res <- make_results(dds_name=dds_name, label=label, type=type) # No contrast when running test == 'LRT'
        check_results(res, lrt_design_data, label, test=test, type=type)
      } else if ((!is.null(test) && test == 'LRT') && (!is.null(type) && type %in% c('apeglm','normal'))) {
        # 'coef' is required for shrinkage type == 'apeglm' and 'apeglm'
        res <- make_results(dds_name=dds_name, label=label, type=type, coef=coef)
        check_results(res, lrt_design_data, label, test=test, type=type)
      } else if (!is.null(test) && type != 'apeglm' && !is.null(type)) {
        res <- make_results(dds_name=dds_name, label=label, type=type, contrast=contrast) # Wald, ashr
        check_results(res, lrt_design_data, label, contrast=contrast, test=test, type=type)
      } else if (is.null(test) && type != 'apeglm' && !is.null(type)) {
        res <- make_results(dds_name=dds_name, label=label, type=type, contrast=contrast)
        check_results(res, lrt_design_data, label, contrast=contrast, test=NULL, type=type)
      } else if (!is.null(test) && is.null(type)) {
        res <- make_results(dds_name=dds_name, label=label, type=NULL, contrast=contrast)
        check_results(res, lrt_design_data, label, contrast=contrast, test=test, type=NULL)
      } else if (is.null(test) && is.null(type)) {
        res <- make_results(dds_name=dds_name, label=label, type=NULL, contrast=contrast)
        check_results(res, lrt_design_data, label, contrast=contrast, test=NULL, type=NULL)
      } else if (!is.null(test) && type == 'apeglm') {
        res <- make_results(dds_name=dds_name, label=label, type=type, coef=coef)
        check_results(res, lrt_design_data, label, coef=coef, test=test, type=type)
      } else if (is.null(test) && type == 'apeglm') {
        res <- make_results(dds_name=dds_name, label=label, type=type, coef=coef)
        check_results(res, lrt_design_data, label, coef=coef, test=NULL, type=type)
      } else {
        stop(paste(label, "was not tested"))
      }
    }) # test_that
  } # for type in shrinkage_types
} # for test in tests

test_that("make_results can handle dds object directly", {
  dds <- dds_list[['dds_wald']]
  # Directly pass the dds object
  results <- make_results(dds_name=dds,
                          label='Direct DDS',
                          type='ashr',
                          contrast=c("condition", "treatment", "control"))

  # Check that the res element is a DESeqResults object
  expect_true(inherits(results$res, "DESeqResults"))
}) # test_that

#test_that("make_results errors when a dds_name is passed and dds_list is missing from .GlobalEnv", {
#  remove(dds_list)
#  expect_error(make_results(dds_name='dds_wald',
#                            label='missing dds_list',
#                            type='ashr',
#                            contrast=c("condition", "treatment", "control")),
#                            "Can't find dds_list in global environment.")
#
#}) # test_that
