#library(DESeq2)
#library(testthat)
#devtools::load_all('../../../../lib/lcdbwf')
config <- lcdbwf:::load_config('../../../../workflows/rnaseq/downstream/config.yaml')
source('test-functions.R')

# --------- Test strip_dotted_version_from_dds() ---------- #
test_that("strip_dotted_version_from_dds works", {
  dds <- DESeq2::makeExampleDESeqDataSet()
  # genes are labeled like "gene1", "gene2", etc; change them to "gene1.1",
  # gene2.2", and so on
  rownames(dds) <- paste(rownames(dds), seq(1000), sep='.')
  expect_error(lcdbwf:::strip_dotted_version_from_dds(dds), "Gene names don't appear to be Ensembl")
  forced <- lcdbwf:::strip_dotted_version_from_dds(dds, force=TRUE)
  expect_equal(rownames(forced)[1], "gene1")
  rownames(dds) <- paste0("ENS", rownames(dds), '.', seq(1000))
  fixed <- lcdbwf:::strip_dotted_version_from_dds(dds)
  expect_equal(rownames(fixed)[1], "ENSgene1")
}) # test_that
# --------------------------------------------------------- #

# --------------------- Test make_dds() ------------------- #
test_that("make_dds handles minimum required design data of sampletale,
           and design formula, using example featurecounts table", {
  design_data <- make_design_data()
  dds <- make_dds(design_data,
                  config=config,
                  featureCounts='featurecounts.txt',
                  parallel=config$parallel$parallel)
  # Check that the dds object is a DESeqDataSet (not NULL)
  expect_true(inherits(dds, "DESeqDataSet"))
  # Check that the WaldStatistic has been generated as Wald is the default test that
  # is run with minimal input to make_dds()
  expect_true(any(grepl('Wald',  colnames(mcols(dds)))))
  expect_false(any(grepl('LRT',  colnames(mcols(dds)))))
}) # test_that

test_that("make_dds handles design data where Wald test is specified explicitly", {
  design_data <- make_design_data()
  design_data$test <- 'Wald'
  dds <- make_dds(design_data,
                  config=config,
                  featureCounts='featurecounts.txt',
                  parallel=config$parallel$parallel)
  # Check that the dds object is a DESeqDataSet (not NULL)
  expect_true(inherits(dds, "DESeqDataSet"))
  # Check that the WaldStatistic has been generated when Wald is passed to make_dds()
  expect_true(any(grepl('Wald',  colnames(mcols(dds)))))
  expect_false(any(grepl('LRT',  colnames(mcols(dds)))))
}) # test_that

test_that("make_dds handles required design data for LRT", {
  design_data <- make_design_data()
  design_data$test <- 'LRT'
  design_data$reduced_design <- ~1
  dds <- make_dds(design_data,
                  config=config,
                  featureCounts='featurecounts.txt',
                  parallel=config$parallel$parallel)
  # Check that the dds object is a DESeqDataSet (not NULL)
  expect_true(inherits(dds, "DESeqDataSet"))
  # Check that the WaldStatistic has been generated
  expect_true(any(grepl('LRT',  colnames(mcols(dds)))))
  expect_false(any(grepl('Wald',  colnames(mcols(dds)))))
}) # test_that

# Now we intentionally call make_results with incompatible parameters
test_that("make_dds errors on invalid 'test' option", {
  design_data <- make_design_data()
  design_data$test <- "invalid_test_option"
  expect_error(make_dds(design_data,
                        config=config,
                        featureCounts='featurecounts.txt',
                        parallel=config$parallel$parallel),
               paste("Valid options for test are \\'Wald\\' \\(default\\) or \\'LRT\\'. You chose,", design_data$test))
}) # test_that

test_that("make_dds errors on missing reduced design when 'test' is set to LRT", {
  design_data <- make_design_data()
  design_data$test <- "LRT"
  expect_error(make_dds(design_data,
                        config=config,
                        featureCounts='featurecounts.txt',
                        parallel=config$parallel$parallel),
               "When using LRT, reduced_design must be provided")
}) # test_that

test_that("make_dds errors on missing test argument when reduced design is provided", {
  design_data <- make_design_data()
  design_data$reduced_design <- ~1
  expect_error(make_dds(design_data,
                        config=config,
                        featureCounts='featurecounts.txt',
                        parallel=config$parallel$parallel),
               "You included a reduced design formula but did not specify test = 'LRT'")
}) # test_that
# --------------------------------------------------- #

# -------------- collapseReplicates2 ---------------- #
test_that("collapseReplicates2 collapses the two control replicates and two treatment replicates to
           a single control row and a single treatment row in colData. Row names should still match
           column 1.", {
  # Setup a DESeqDataSet with replicates
  design_data <- make_design_data()
  # Add bio_rep to colData
  dds <- make_dds(design_data,
                  config=config,
                  featureCounts='featurecounts.txt',
                  parallel=config$parallel$parallel)

  # Perform replicate collapsing
  dds_collapsed <- collapseReplicates2(dds, dds$group)
  colData_collapsed <- as.data.frame(colData(dds_collapsed))
  expect_equal(length(as.character(colData_collapsed$group[colData_collapsed$group == 'control'])), 1)
  expect_equal(length(as.character(colData_collapsed$group[colData_collapsed$group == 'treatment'])), 1)
  # Check if the first column of colData matches rownames
  expect_equal(rownames(colData_collapsed), colData_collapsed[,1])
}) # test_that
# --------------------------------------------------- #
