library(DESeq2)
library(testthat)
library(rlang)
library(stringr)
library(BiocParallel)
library(future)
devtools::load_all('../../../../lib/lcdbwf')
config <- lcdbwf:::load_config('config.yaml')
source('test-functions.R')
register(MulticoreParam(workers = future::availableCores()), default=TRUE)
#param <- bpparam()
#number_of_cores <- param$workers
#print(number_of_cores)

# Test all combinations of test and type
# NULL shrinkage type skips lfcShrink
# NULL test type runs Wald (default test)
tests <- list('Wald', 'LRT', NULL)
shrinkage_types <- list('ashr', 'apeglm', 'normal', NULL)
contrast <- c("condition", "treatment", "control")
coef <- "condition_treatment_vs_control"
# Make the dds_list containing dds_wald and dds_lrt dds objects
# Also save the full and reduced design formulas used to create dds_lrt
dds_and_lrt_design <- make_lists()
dds_list <- dds_and_lrt_design$dds_list # The get_dds call in make_results requires dds_list to be in .GlobalEnv
lrt_design_data <- dds_and_lrt_design$lrt_design_data
#test_make_results(tests, shrinkage_types, contrast, coef, dds_list, lrt_design_data)

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
        res <- make_results(dds_name=dds_name, label=label, test=test, type=NULL) # No contrast when running test == 'LRT'
        check_results(res, lrt_design_data, label, test=test, type=NULL)
      } else if ((!is.null(test) && test == 'LRT') && (!is.null(type) && !type %in% c('apeglm','normal'))) {
        res <- make_results(dds_name=dds_name, label=label, test=test, type=type) # No contrast when running test == 'LRT'
        check_results(res, lrt_design_data, label, test=test, type=type)
      } else if ((!is.null(test) && test == 'LRT') && (!is.null(type) && type %in% c('apeglm','normal'))) {
        # No contrast when running test == 'LRT'. But coef is required for shrinkage type == 'apeglm' and 'apeglm'
        res <- make_results(dds_name=dds_name, label=label, test=test, type=type, coef=coef)
        check_results(res, lrt_design_data, label, test=test, type=type)
      } else if (!is.null(test) && type != 'apeglm' && !is.null(type)) {
        res <- make_results(dds_name=dds_name, label=label, test=test, type=type, contrast=contrast) # Wald, ashr
#print("str(res) ---------------- ")
#print(str(res))
#print("str(res$res) ---------------- ")
#print(str(res$res))
#print("names of metadata of res: -----------")
#print(names(metadata(res$res)))
#print("$type of metadata of res: ----------")
#print(metadata(res$res)$type)
#print("Expected type: ------------")
#print(type)
        check_results(res, lrt_design_data, label, contrast=contrast, test=test, type=type)
      } else if (is.null(test) && type != 'apeglm' && !is.null(type)) {
        res <- make_results(dds_name=dds_name, label=label, test=NULL, type=type, contrast=contrast)
        check_results(res, lrt_design_data, label, contrast=contrast, test=NULL, type=type)
      } else if (!is.null(test) && is.null(type)) {
        res <- make_results(dds_name=dds_name, label=label, test=test, type=NULL, contrast=contrast)
        check_results(res, lrt_design_data, label, contrast=contrast, test=test, type=NULL)
      } else if (is.null(test) && is.null(type)) {
        res <- make_results(dds_name=dds_name, label=label, test=NULL, type=NULL, contrast=contrast)
        check_results(res, lrt_design_data, label, contrast=contrast, test=NULL, type=NULL)
      } else if (!is.null(test) && type == 'apeglm') {
        res <- make_results(dds_name=dds_name, label=label, test=test, type=type, coef=coef)
        check_results(res, lrt_design_data, label, coef=coef, test=test, type=type)
      } else if (is.null(test) && type == 'apeglm') {
        res <- make_results(dds_name=dds_name, label=label, test=NULL, type=type, coef=coef)
        check_results(res, lrt_design_data, label, coef=coef, test=NULL, type=type)
      } else {
        stop(paste(label, "was not tested"))
      }
    }) # test_that make_results works correctly with each combination of test and type
  } # for type in shrinkage_types
} # for test in tests

# Now we intentionally call make_results with incompatible parameters
test_that("make_results errors on invalid 'test' option", {
  design_data <- make_design_data()
  design_data$test <- "invalid_test_option"

        res <- make_results(dds_name=dds_name, label=label, test=NULL, type=NULL, contrast=contrast)
  expect_error(make_dds(design_data,
                        config=config,
                        featureCounts='featurecounts.txt',
                        parallel=config$parallel$parallel),
              paste("Valid options for test are \\'Wald\\' \\(default\\) or \\'LRT\\'. You chose,", test))
}) # test_that make_dds errors on invalid test option

test_that("make_results can handle dds object directly", {
  design_data <- make_design_data()
  design_data$test <- 'Wald'
  make_featurecounts_file()
  dds <- make_dds(design_data,
                  config=config,
                  featureCounts='featurecounts.txt',
                  parallel=config$parallel$parallel)

  # Directly pass the dds object
  results <- make_results(dds_name=dds,
                          label='Direct DDS',
                          test='Wald',
                          type='ashr',
                          contrast=c("condition", "treated", "control"))

  # Check that the res element is a DESeqResults object
  expect_true(inherits(results$res, "DESeqResults"))
  # Check that the metadata of the results object includes the correct type
  expect_true(metadata(results$res)$type == "ashr")
})

test_that("make_results handles missing 'samplename' column", {
  design_data <- make_design_data()
  design_data$test <- 'Wald'
  make_featurecounts_file()
  dds <- make_dds(design_data,
                  config=config,
                  featureCounts='featurecounts.txt',
                  parallel=config$parallel$parallel)

  # Remove the 'samplename' column to trigger error
  colData(dds)$samplename <- NULL

  expect_error(dds_coefs(dds, colour=='white'),
               "Need to have 'samplename' as a column in colData")
})

