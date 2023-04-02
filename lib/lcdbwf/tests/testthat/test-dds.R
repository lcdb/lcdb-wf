test_that("stripping gene versions", {
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

})
