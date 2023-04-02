sampletable <- data.frame(
        samplename=c('sample1', 'sample2', 'sample3', 'sample4'),
        group=factor(c('control', 'control', 'treatment', 'treatment'))
        )


test_that("load from combined featurecounts", {
    dds <- lcdbwf:::DESeqDataSetFromCombinedFeatureCounts(
        "featurecounts.txt",
        sampletable,
        design=~group
    )
    expect_s4_class(dds, "DESeqDataSet")
})

test_that("load from combined featurecounts, subset sampletable, no subset.counts", {
  expect_error(
    lcdbwf:::DESeqDataSetFromCombinedFeatureCounts(
      "featurecounts.txt",
      sampletable %>% dplyr::filter(samplename != "sample4"),
      design=~group),
    paste(
      "The following samples are in the counts data but not the sampletable.",
      "If this is intended, consider using `subset.counts=TRUE` to remove",
      "them from the counts: sample4")
    )
})

test_that("load from combined featurecounts, subset sampletable, with subset.counts", {
    dds <- lcdbwf:::DESeqDataSetFromCombinedFeatureCounts(
        "featurecounts.txt",
        sampletable %>% dplyr::filter(samplename != "sample4"),
        design=~group, subset_counts=TRUE)
    expect_s4_class(dds, "DESeqDataSet")
})


test_that("load from combined featurecounts, using custom func", {
  sampletable2 <- sampletable %>%
    dplyr::mutate(samplename=stringr::str_replace(samplename, "sample", "smp"))

  custom.func <- function(x) {
    x <- x %>%
        stringr::str_remove_all('data/rnaseq_samples/') %>%
        stringr::str_remove_all('.cutadapt.bam') %>%
        stringr::str_replace("sample", "smp") %>%
        stringr::str_split(stringr::fixed('/'), simplify=TRUE)
    x[,1]
}
  dds <- lcdbwf:::DESeqDataSetFromCombinedFeatureCounts(
      "featurecounts.txt",
      sampletable2,
      sample_func=custom.func,
      design=~group)

    expect_s4_class(dds, "DESeqDataSet")
})
