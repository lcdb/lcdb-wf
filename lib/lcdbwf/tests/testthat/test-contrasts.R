devtools::load_all('../../../../lib/lcdbwf')

# Used for the %||% operator
library(rlang)

#register(MulticoreParam(config$parallel$cores))
config <- lcdbwf:::load_config('../../../../workflows/rnaseq/downstream/config.yaml')
source('test-functions.R')

is_deseq_res <- function(x) {
  inherits(x, "DESeqResults")
}
# Test all combinations of test and type
# NULL shrinkage type skips lfcShrink
# NULL test type runs Wald (default test)
tests <- list('Wald', 'LRT', NULL)
shrinkage_types <- list('ashr', 'apeglm', 'normal', NULL)
contrast <- c("group", "treatment", "control")
coef <- "group_treatment_vs_control"
dds_list <- make_dds_list(config)

# Ensure dds_list makes it into the global environment, no matter what fancy
# stuff {testthat} is doing.
assign("dds_list", dds_list, envir=.GlobalEnv)

lrt_design_data <- make_lrt_design_data()

# Each row in the ASCII table indicates which combination of test, type, coef, and contrast
# is tested by the respective indexed conditional statement in the following test_that code.

#+---------+-------+-------+------+----------+-------+
#| Results | Test  | Type  | Coef | Contrast | Check |
#+---------+-------+-------+------+----------+-------+
#|    1    | LRT   | NULL  | -    | -        |   E   |
#|    2    | LRT   | ashr  | -    | -        |   F   |
#|    2    | LRT   | apeglm| -    | -        |   F   |
#|    2    | LRT   | normal| -    | -        |   F   |
#|    5    | Wald  | NULL  | -    | yes      |   C   |
#|    3    | Wald  | ashr  | -    | yes      |   A   |
#|    7    | Wald  | apeglm| yes  | -        |   B   |
#|    3    | Wald  | normal| -    | yes      |   C   |
#|    6    | NULL  | NULL  | -    | yes      |   C   |
#|    4    | NULL  | ashr  | -    | yes      |   A   |
#|    8    | NULL  | apeglm| yes  | -        |   B   |
#|    4    | NULL  | normal| -    | yes      |   C   |
#+---------+-------+-------+------+----------+-------+

for (test in tests) {
  for (type in shrinkage_types) {
    if (test == 'Wald' || is.null(test)) {
      dds_name <- 'dds_wald'
    } else if (test == 'LRT') {
      dds_name <- 'dds_lrt'
    }
    label <- paste0("test=", test %||% "NULL/default (Wald)", ", type=", type %||% "NULL (Skip)")
    test_that(paste("make_results works correctly with type =", type, "and", test, "being detected automatically
                     from DDS"), {
      # 'Results' from the table above
      # 1
      if ((!is.null(test) && test == 'LRT') && is.null(type)) {
        res <- make_results(dds_name=dds_name, label=label, type=NULL)
      # 2
      } else if ((!is.null(test) && test == 'LRT') && (!is.null(type) && type %in% c('ashr', 'apeglm', 'normal'))) {
      # 'Check' from the table above
      # F
        expect_error(make_results(dds_name=dds_name, label=label, type=type),
          "You cannot pass a non-NULL or missing type to make_results with an LRT dds object. For LRT, LFC values are set to 0 and should not be passed to lfcShrink. Use type == NULL in make_results for LRT DDS objects.")
        return()
      # 3
      } else if (!is.null(test) && type != 'apeglm' && !is.null(type)) {
        res <- make_results(dds_name=dds_name, label=label, type=type, contrast=contrast)
      # 4
      } else if (is.null(test) && type != 'apeglm' && !is.null(type)) {
        res <- make_results(dds_name=dds_name, label=label, type=type, contrast=contrast)
      # 5
      } else if (!is.null(test) && is.null(type)) {
        res <- make_results(dds_name=dds_name, label=label, type=NULL, contrast=contrast)
      # 6
      } else if (is.null(test) && is.null(type)) {
        res <- make_results(dds_name=dds_name, label=label, type=NULL, contrast=contrast)
      # 7
      } else if (!is.null(test) && type == 'apeglm') {
        res <- make_results(dds_name=dds_name, label=label, type=type, coef=coef)
      # 8
      } else if (is.null(test) && type == 'apeglm') {
        res <- make_results(dds_name=dds_name, label=label, type=type, coef=coef)
      } else {
        stop(paste(label, "was not tested"))
      }

      # Check make_results output for each possible combination of test and type
      expect_true(is_deseq_res(res$res))
      expect_true(identical(names(res), c('res', 'dds', 'label')))
      lrt_mcols_description <- paste0(as.character(lrt_design_data$design)[1], " ",
                                      as.character(lrt_design_data$design)[2], "' vs '",
                                      as.character(lrt_design_data$reduced_design)[1], " ",
                                      as.character(lrt_design_data$reduced_design)[2], "'")
      # 'Check' from the table above
      # A
      if ((is.null(test) || test == 'Wald') && (!is.null(type) && type == 'ashr')) {
        expected_char <- paste(test %||% 'Wald', "test p-value:", contrast[1], contrast[2], "vs", contrast[3])
        expect_true(S4Vectors::mcols(res$res)$description[4] == expected_char)
      # B
      } else if ((is.null(test) || test == 'Wald') && (!is.null(type) && type == 'apeglm')) {
        coef <- stringr::str_split(coef, "_")[[1]]
        expected_char <- paste(test %||% 'Wald', "test p-value:", coef[1], coef[2], coef[3], coef[4])
        expect_true(S4Vectors::mcols(res$res)$description[4] == expected_char)
      # C
      } else if ((is.null(test) || test == 'Wald') && (is.null(type) || type == 'normal')) {
        expected_char <- paste(test %||% 'Wald', "statistic:", contrast[1], contrast[2], "vs", contrast[3])
        expect_true(S4Vectors::mcols(res$res)$description[4] == expected_char)
      # D
      } else if ((!is.null(test) && test == 'LRT') && (!is.null(type) && type != 'normal')) {
        expected_char <- paste0(test, " p-value: '", lrt_mcols_description)
        expect_true(S4Vectors::mcols(res$res)$description[4] == expected_char)
      # E
      } else if ((!is.null(test) && test == 'LRT') && (is.null(type) || type == 'normal')) {
        expect_true(all(res$res$log2FoldChange == 0))
        expected_char <- paste0(test, " statistic: '", lrt_mcols_description)
        expect_true(S4Vectors::mcols(res$res)$description[4] == expected_char)
      } else {
        stop(paste(label, 'was not checked'))
      }
      # Check for expected type stored in the result's metadata
      if (!is.null(type)) {
        expect_true(identical(S4Vectors::metadata(res$res)$type, type))
      } else if (is.null(type)) {
        expect_true(is.null(S4Vectors::metadata(res$res)$type))
      }
    }) # test_that
  } # for type in shrinkage_types
} # for test in tests


# -------------------------- direct dds -------------------------- #
test_that("make_results can handle dds object directly", {
  dds <- dds_list[['dds_wald']]
  # Directly pass the dds object
  res <- make_results(dds_name=dds,
                      label='Direct DDS',
                      type='ashr',
                      contrast=contrast)
  expect_true(is_deseq_res(res$res))
  expect_true(identical(names(res), c('res', 'dds', 'label')))
  expect_true(any(grepl('Wald', S4Vectors::mcols(res$res)$description[4])))
}) # test_that
# ---------------------------------------------------------------- #


# ---- make_results with dds_lrt but with wald test specified ---- #
# Similar structure to the ASCII table depicted tests from above
# but with the 'test' argument included in make_results
# and test == 'LRT'
test_that("make_results errors when user passes mismatched test == 'Wald' with LRT DDS", {
  expect_error(make_results(dds_name='dds_lrt',
                            label='Shrink lrt results',
                            type=NULL,
                            test='Wald'),
                            "The 'test' passed to make_results was set to 'Wald' but 'LRT' has been detected in dds")
}) # test_that
# ---------------------------------------------------------------- #


# ---- make_results with dds_wald but with LRT test specified ---- #
# Similar structure to the ASCII table depicted tests from above
# but with the 'test' argument included in make_results
# and test == 'LRT'
test_that("make_results errors when user passes mismatched test == 'LRT' with Wald DDS", {
  expect_error(make_results(dds_name='dds_wald',
                            label='Shrink lrt results',
                            type=NULL,
                            test='LRT'),
                            "The 'test' passed to make_results was set to 'LRT' but 'Wald' has been detected in dds")
}) # test_that
# ---------------------------------------------------------------- #


# ----- Attempt to shrink LRT results but with test included ----- #
# Similar structure to the ASCII table depicted tests from above
# but with the 'test' argument included in make_results
# and test == 'LRT'
for (type in c('ashr', 'apeglm', 'normal')) {
  test_that("make_results errors when user attempts to run lfcShrink by defining a non-NULL or missing type when test == 'LRT'", {
    expect_error(make_results(dds_name='dds_lrt',
                              label='Shrink lrt results',
                              type=type,
                              test='LRT'),
                              "You cannot pass a non-NULL or missing type to make_results with test == 'LRT'. For LRT, LFC values are set to 0 and should not be passed to lfcShrink. Use type == NULL in make_results for LRT DDS objects.")
  }) # test_that
} # for type
# ---------------------------------------------------------------- #


# -------------- missing both test and type with LRT ------------- #
# Similar structure to the ASCII table depicted tests from above
# with the 'test' argument also missing in make_results
# and test also missing. When test is missing, test should be detected from dds_lrt as 'LRT'.
# With type also missing, type should be set as the current default: 'ashr' as of 5-30-2024.
# This combination of test and type is incompatible and so the following error message should
# be returned.
test_that("make_results errors when user attempts to run lfcShrink by defining a missing type when test == 'LRT'", {
  expect_error(make_results(dds_name='dds_lrt',
                            label='missing test and type of LRT DDS'),
                            "You cannot pass a non-NULL or missing type to make_results with an LRT dds object. For LRT, LFC values are set to 0 and should not be passed to lfcShrink. Use type == NULL in make_results for LRT DDS objects.")
}) # test_that
# ---------------------------------------------------------------- #


# ------------- test included for all types dds Wald ------------- #
# Similar structure to the ASCII table depicted tests from above
# but with the 'test' argument included in make_results
# and test == 'Wald'
for (type in c('ashr', 'apeglm', 'normal')) {
  test_that("make_results returns a DESeqResults object with !all res$res$LFC == 0 when user passes a defined type along with test == 'Wald'", {
    if (type != 'apeglm') {
      res <- make_results(dds_name='dds_wald',
                          label='Shrink Wald results',
                          contrast=contrast,
                          type=type,
                          test='Wald')
      expect_true(is_deseq_res(res$res))
      expect_true(identical(names(res), c('res', 'dds', 'label')))
      expect_true(!all(res$res$log2FoldChange == 0))
      expect_true(any(grepl('Wald', S4Vectors::mcols(res$res)$description[4])))
    } else if (type == 'apeglm') {
      res <- make_results(dds_name='dds_wald',
                                        label='Shrink Wald results',
                                        coef=coef,
                                        type=type,
                                        test='Wald')
      expect_true(is_deseq_res(res$res))
      expect_true(identical(names(res), c('res', 'dds', 'label')))
      expect_true(!all(res$res$log2FoldChange == 0))
      expect_true(any(grepl('Wald', S4Vectors::mcols(res$res)$description[4])))
    } # if type
  }) # test_that
} # for type
# ---------------------------------------------------------------- #


# ---------------------- missing dds_list ------------------------ #
orig_dds_list <- dds_list
remove(dds_list, envir=.GlobalEnv)
test_that("make_results errors when a dds_name is passed and dds_list is missing from .GlobalEnv", {
  expect_error(make_results(dds_name='dds_wald',
                            label='missing dds_list',
                            type='ashr',
                            contrast=contrast),
                            "Can't find dds_list in global environment.")
}) # test_that
# Put it back into the global env
assign("dds_list", orig_dds_list, envir=.GlobalEnv)
# ---------------------------------------------------------------- #
