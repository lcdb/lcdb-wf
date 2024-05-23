library(DESeq2)
library(testthat)
library(rlang)
library(stringr)
devtools::load_all('../../../../lib/lcdbwf')
config <- lcdbwf:::load_config('config.yaml')
source('test-functions.R')

# Test all combinations of test and type
# NULL shrinkage type skips lfcShrink
# NULL test type runs Wald (default test)
tests <- list('Wald', 'LRT', NULL)
shrinkage_types <- list('ashr', 'apeglm', 'normal', NULL)
contrast <- c("condition", "treatment", "control")
coef <- "condition_treatment_vs_control"
# Make the dds_list containing dds_wald and dds_lrt dds objects
# Also save the full and reduced design formulas used to create dds_lrt
#dds_and_lrt_design <- make_lists()
#dds_list <- dds_and_lrt_design$dds_list # The get_dds call in make_results requires dds_list to be in .GlobalEnv
#lrt_design_data <- dds_and_lrt_design$lrt_design_data
test_make_results(tests, shrinkage_types, contrast, coef, dds_list, lrt_design_data)

# Now we intentionally call make_results with incompatible sets of parameters
# based on what I think is likely
test_that("make_results errors on invalid 'test' option", {
  design_data <- make_design_data()
  design_data$test <- "invalid_test_option"
  design_data$reduced_design <- ~1

  expect_error(make_dds(design_data,
                        config=config,
                        featureCounts='featurecounts.txt',
                        parallel=config$parallel$parallel),
               "Valid options for test are \\'Wald\\' \\(default\\) or \\'LRT\\'")
})

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

