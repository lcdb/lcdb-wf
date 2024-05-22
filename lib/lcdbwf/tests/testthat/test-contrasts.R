library(DESeq2)
devtools::load_all('../../../../lib/lcdbwf')
config <- lcdbwf:::load_config('config.yaml')
library(testthat)
source('test-functions.R')
library(rlang)
library(stringr)

# Helper function to run make_results and check the output
make_results_and_check <- function(dds_name=dds_name,
                                   label=label,
                                   test=NULL,
                                   type=NULL,
                                   lrt_design_data.=lrt_design_data,
                                   contrast=NULL,
                                   coef=NULL) {
  print(label)
  if (type != 'apeglm') {
    # Use the 'contrast' argument when shrinkage type != 'apeglm'
    res <- make_results(dds_name=dds_name,
                        label=label,
                        test=test,
                        type=type,
                        contrast=contrast)
  } else if (type == 'apeglm') {
    # Use the 'coef' argument instead of 'contrast' when shrinkage type == 'apeglm'
    res <- make_results(dds_name=dds_name,
                        label=label,
                        test=test,
                        type=type,
                        coef=coef)
  }

  # Check that res were returned by make_res
  expect_true(!is.null(res))
  # Check that the res element returned by make_res is a DESeqres object
  expect_true(inherits(res$res, "DESeqResults"))
  # Save a character representing the source of LRT pvalue for comparison with make_results output
  lrt_mcols_description <- paste0(as.character(lrt_design_data.$design)[1], " ",
                                  as.character(lrt_design_data.$design)[2], "' vs '",
                                  as.character(lrt_design_data.$reduced_design)[1], " ",
                                  as.character(lrt_design_data.$reduced_design)[2], "'")

  # Check the metadata in res for correct test and coef/contrast
  # based on each combination of test, type and coef/contrast arguments
  if ((is.null(test) || test == 'Wald') && type == 'apeglm') {
    coef <- str_split(coef, "_")[[1]]
    expected_char <- paste(test, "test p-value:", coef[1], coef[2], coef[3], coef[4])
    expect_true(mcols(res$res)$description[4] == expected_char)
  # Check that res for the correct test was extracted for all tests and types excluding 'apeglm'
  } else if ((is.null(test) || test == 'Wald') && type != 'normal') {
    expected_char <- paste(test, "test p-value:", contrast[1], contrast[2], "vs", contrast[3])
    expect_true(mcols(res$res)$description[4] == expected_char)
  } else if ((is.null(test) || test == 'Wald')  && type == 'normal') {
    expected_char <- paste(test, "statistic:", contrast[1], contrast[2], "vs", contrast[3])
    expect_true(mcols(res$res)$description[4] == expected_char)
  } else if (test == 'LRT' && type != 'normal') {
    expected_char <- paste0(test, " p-value: '", lrt_mcols_description)
    expect_true(mcols(res$res)$description[4] == expected_char)
  } else if (test == 'LRT' && type == 'normal') {
    expected_char <- paste0(test, " statistic: '", lrt_mcols_description)
    expect_true(mcols(res$res)$description[4] == expected_char)
  } else {
    stop(paste(label, 'was not checked'))
  }

  #   mdcat(mcols(res_i)$description[2]) # For log fold change value possibly want to test for LRT
  # Check that the metadata of the res object includes the correct shrinkage type
  if (!is.null(type)) {
    expect_true(metadata(res$res)$type == type)
  } else if (is.null(type)) {
    expect_true(is.null(metadata(res$res)$type))
  }
} # make_results_and_check

# Get the dds_list containing dds_wald and dds_lrt dds objects
# Also save the full and reduced design formulas used to create dds_lrt
dds_and_lrt_design <- make_dds_list()
dds_list <- dds_and_lrt_design$dds_list
lrt_design_data <- dds_and_lrt_design$lrt_design_data

# Test all combinations of test and type
# NULL shrinkage type skips lfcshrink
# NULL test type runs Wald test (default test)
tests <- list('Wald', 'LRT', NULL)
shrinkage_types <- list(NULL, 'ashr', 'apeglm', 'normal')

for (test in tests) {
  for (type in shrinkage_types) {
    if (test == 'Wald' || is.null(test)) {
      dds_name <- 'dds_wald'
    } else if (test == 'LRT') {
      dds_name <- 'dds_lrt'
    }
    test_label <- paste0("test=", ifelse(is.null(test), "NULL/default (Wald)", test), ", type=", ifelse(is.null(type), "NULL (Skip)", type))
    test_that(paste("make_results works correctly with", test_label), {
      if (type != 'apeglm') {
        make_results_and_check(dds_name=dds_name,
                               label=test_label,
                               test=test,
                               type=type,
                               lrt_design_data=lrt_design_data,
                               contrast=c("condition", "treatment", "control"))
      } else if (type == 'apeglm') {
        make_results_and_check(dds_name=dds_name,
                               label=test_label,
                               test=test,
                               type=type,
                               lrt_design_data,
                               coef="condition_treatment_vs_control")
      } # else if type == 'apeglm'
    }) # test_that make_results works correctly with each combination of test and type
  } # for type in shrinkage_types
} # for test in tests

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

