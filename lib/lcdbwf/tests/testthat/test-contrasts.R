library(DESeq2)
library(testthat)
library(rlang)
library(stringr)
library(BiocParallel)
devtools::load_all('../../../../lib/lcdbwf')
config <- lcdbwf:::load_config('../../../../workflows/rnaseq/downstream/config.yaml')
source('test-functions.R')
register(MulticoreParam(config$parallel$cores))

# Function to create design data for LRT test
make_lrt_design_data <- function() {
  lrt_design_data <- make_design_data()
  lrt_design_data$test <- 'LRT'
  lrt_design_data$reduced_design <- ~1
  return(lrt_design_data)
} # make_lrt_design_data

make_dds_list <- function() {
  # Create design data and dds object for Wald test type
  wald_design_data <- make_design_data()
  dds_wald <- make_dds(wald_design_data,
                       config=config,
                       featureCounts='featurecounts.txt',
                       parallel=config$parallel$parallel)

  lrt_design_data <- make_lrt_design_data()
  dds_lrt <- make_dds(lrt_design_data,
                      config=config,
                      featureCounts='featurecounts.txt',
                      parallel=config$parallel$parallel)

  # Create dds_list
  dds_list <- list(dds_wald=dds_wald, dds_lrt=dds_lrt)
  return(dds_list)
} # make_dds_list

# Test all combinations of test and type
# NULL shrinkage type skips lfcShrink
# NULL test type runs Wald (default test)
tests <- list('Wald', 'LRT', NULL)
shrinkage_types <- list('ashr', 'apeglm', 'normal', NULL)
contrast <- c("group", "treatment", "control")
coef <- "group_treatment_vs_control"
dds_list <- make_dds_list()
lrt_design_data <- make_lrt_design_data()

### TESTING ###
#test <- 'Wald'
#type <- 'ashr'
#dds_name <- 'dds_wald'
#contrast <- c("condition", "treatment", "control")
#label <- paste0("test=", test %||% "NULL/default (Wald)", ", type=", type %||% "NULL (Skip)")
##############

# Each row in the ASCII table indicates which combination of test, type, coef, and contrast
# is tested by the respective indexed conditional statement in the following test_that code.

#+---------+-------+-------+------+----------+-------+
#| Results | Test  | Type  | Coef | Contrast | Check |
#+---------+-------+-------+------+----------+-------+
#|    1    | LRT   | NULL  | -    | -        |   E   |
#|    2    | LRT   | ashr  | -    | -        |   D   |
#|    3    | LRT   | apeglm| yes  | -        |   D   |
#|    3    | LRT   | normal| yes  | -        |   E   |
#|    6    | Wald  | NULL  | -    | yes      |   C   |
#|    4    | Wald  | ashr  | -    | yes      |   A   |
#|    8    | Wald  | apeglm| yes  | -        |   B   |
#|    4    | Wald  | normal| -    | yes      |   C   |
#|    7    | NULL  | NULL  | -    | yes      |   C   |
#|    5    | NULL  | ashr  | -    | yes      |   A   |
#|    9    | NULL  | apeglm| yes  | -        |   B   |
#|    5    | NULL  | normal| -    | yes      |   C   |
#+---------+-------+-------+------+----------+-------+

for (test in tests) {
  for (type in shrinkage_types) {
    if (test == 'Wald' || is.null(test)) {
      dds_name <- 'dds_wald'
    } else if (test == 'LRT') {
      dds_name <- 'dds_lrt'
    }
    label <- paste0("test=", test %||% "NULL/default (Wald)", ", type=", type %||% "NULL (Skip)")
    test_that(paste("make_results works correctly with", label), {
      # 'Results' from the table above
      # 1
      if ((!is.null(test) && test == 'LRT') && is.null(type)) {
        res <- make_results(dds_name=dds_name, label=label, type=NULL)
      # 2
      } else if ((!is.null(test) && test == 'LRT') && (!is.null(type) && !type %in% c('apeglm','normal'))) {
        res <- make_results(dds_name=dds_name, label=label, type=type)
      # 3
      } else if ((!is.null(test) && test == 'LRT') && (!is.null(type) && type %in% c('apeglm','normal'))) {
        res <- make_results(dds_name=dds_name, label=label, type=type, coef=coef)
      # 4
      } else if (!is.null(test) && type != 'apeglm' && !is.null(type)) {
        res <- make_results(dds_name=dds_name, label=label, type=type, contrast=contrast)
      # 5
      } else if (is.null(test) && type != 'apeglm' && !is.null(type)) {
        res <- make_results(dds_name=dds_name, label=label, type=type, contrast=contrast)
      # 6
      } else if (!is.null(test) && is.null(type)) {
        res <- make_results(dds_name=dds_name, label=label, type=NULL, contrast=contrast)
      # 7
      } else if (is.null(test) && is.null(type)) {
        res <- make_results(dds_name=dds_name, label=label, type=NULL, contrast=contrast)
      # 8
      } else if (!is.null(test) && type == 'apeglm') {
        res <- make_results(dds_name=dds_name, label=label, type=type, coef=coef)
      # 9
      } else if (is.null(test) && type == 'apeglm') {
        res <- make_results(dds_name=dds_name, label=label, type=type, coef=coef)
      } else {
        stop(paste(label, "was not tested"))
      }

      # Check make_results output for each possible combination of test and type
      print(label)
      expect_true(inherits(res$res, "DESeqResults"))
      expect_true(identical(names(res), c('res', 'dds', 'label')))
      lrt_mcols_description <- paste0(as.character(lrt_design_data$design)[1], " ",
                                      as.character(lrt_design_data$design)[2], "' vs '",
                                      as.character(lrt_design_data$reduced_design)[1], " ",
                                      as.character(lrt_design_data$reduced_design)[2], "'")
      # 'Check' from the table above
      # A
      if ((is.null(test) || test == 'Wald') && (!is.null(type) && type == 'ashr')) {
        expected_char <- paste(test %||% 'Wald', "test p-value:", contrast[1], contrast[2], "vs", contrast[3])
        expect_true(mcols(res$res)$description[4] == expected_char)
      # B
      } else if ((is.null(test) || test == 'Wald') && (!is.null(type) && type == 'apeglm')) {
        coef <- str_split(coef, "_")[[1]]
        expected_char <- paste(test %||% 'Wald', "test p-value:", coef[1], coef[2], coef[3], coef[4])
        expect_true(mcols(res$res)$description[4] == expected_char)
      # C
      } else if ((is.null(test) || test == 'Wald') && (is.null(type) || type == 'normal')) {
        expected_char <- paste(test %||% 'Wald', "statistic:", contrast[1], contrast[2], "vs", contrast[3])
        expect_true(mcols(res$res)$description[4] == expected_char)
      # D
      } else if ((!is.null(test) && test == 'LRT') && (!is.null(type) && type != 'normal')) {
        expected_char <- paste0(test, " p-value: '", lrt_mcols_description)
        expect_true(mcols(res$res)$description[4] == expected_char)
      # E
      } else if ((!is.null(test) && test == 'LRT') && (is.null(type) || type == 'normal')) {
        expected_char <- paste0(test, " statistic: '", lrt_mcols_description)
        expect_true(mcols(res$res)$description[4] == expected_char)
      } else {
        stop(paste(label, 'was not checked'))
      }
      # Check for expected type stored in the result's metadata
      if (!is.null(type)) {
        expect_true(identical(metadata(res$res)$type, type))
      } else if (is.null(type)) {
        expect_true(is.null(metadata(res$res)$type))
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
                          contrast=c("group", "treatment", "control"))

  # Check that the res element is a DESeqResults object
  expect_true(inherits(results$res, "DESeqResults"))
}) # test_that

test_that("make_results errors when user attempts to run lfcShrink by defining a non-NULL type when test == 'LRT'", {
  expect_error(make_results(dds_name='dds_lrt',
                            label='Shrink lrt results',
                            type='ashr',
                            test == 'LRT'),
                            "You cannot pass a non-NULL type to make_results with test == 'LRT'.
                            For LRT, LFC values are set to 0 and should not be passed to lfcShrink.
                            Use type == NULL in make_results for LRT DDS objects.")
}) # test_that

test_that("make_results errors when user attempts to run lfcShrink by defining a non-NULL type when test is missing", {
  expect_error(make_results(dds_name='dds_lrt',
                            label='Shrink lrt results',
                            type='ashr'),
                            "You cannot pass a non-NULL type to make_results with an LRT dds object.
                            For LRT, LFC values are set to 0 and should not be passed to lfcShrink.
                            Use type == NULL in make_results for LRT DDS objects.")
}) # test_that

test_that("make_results returns a DESeqResults object with all res$res$LFC == 0 when user
           passes type == 'NULL' along with test == 'LRT'", {
  res <- expect_silent(make_results(dds_name='dds_lrt',
                            label='Shrink lrt results',
                            type=NULL,
                            test='LRT'))
      expect_true(inherits(res$res, "DESeqResults"))
      expect_true(identical(names(res), c('res', 'dds', 'label')))
      expect_true(all(res$res$log2FoldChange == 0))
}) # test_that

test_that("make_results returns a DESeqResults object with all res$res$LFC == 0 when user
           passes type == 'NULL' along with missing test parameter'", {
  res <- expect_silent(make_results(dds_name='dds_lrt',
                            label='Shrink lrt results',
                            type=NULL))
      expect_true(inherits(res$res, "DESeqResults"))
      expect_true(identical(names(res), c('res', 'dds', 'label')))
      expect_true(all(res$res$log2FoldChange == 0))
}) # test_that

remove(dds_list)
test_that("make_results errors when a dds_name is passed and dds_list is missing from .GlobalEnv", {
  expect_error(make_results(dds_name='dds_wald',
                            label='missing dds_list',
                            type='ashr',
                            contrast=c("group", "treatment", "control")),
                            "Can't find dds_list in global environment.")

}) # test_that
