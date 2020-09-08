## -----------------------------------------------------------------------------
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(DEGreport)
library(ggplot2)
library(ggrepel)
library(heatmaply)
library(readr)
library(stringr)
library(tibble)
library(purrr)
library(tidyr)

#' Get the OrgDb for the specified organism, using the cached AnnotationHub.
#'
#' @param species Case-sensitive genus and species
#' @param cache Directory in which the AnnotationHub cache is stored
#' @param annotation_key_override If not NA, forces the hub to use this
#'        accession. Use this when you know exactly which OrgDb to use.
#'
#' @return OrgDb object
get.orgdb <- function(species, cache, annotation_key_override=NA){

    # Workaround to allow AnnotationHub to use proxy. See
    # https://github.com/Bioconductor/AnnotationHub/issues/4, and thanks
    # Wolfgang!
    proxy <- Sys.getenv('http_proxy')
    if (proxy == ""){
        proxy <- NULL
    }

    ah <- AnnotationHub(hub=getAnnotationHubOption('URL'),
             cache=cache,
             proxy=proxy,
             localHub=FALSE)

    find.annotationhub.name <- function(species.name, override.code) { #autodetect ah names based on loaded database
        if (is.na(override.code)) {
        ah.query <- query(ah, "OrgDb")
        ah.query.speciesmatch <- grepl(paste("^", species.name, "$", sep=""), ah.query$species)
        ah.query.which <- which(ah.query.speciesmatch)
        stopifnot(length(ah.query.which) > 0) #require at least one match
        if (length(ah.query.which) > 1) { #warn of duplicate matches
           print("WARNING: found multiple candidate species in AnnotationHub: ");
           print(ah.query.speciesmatch)
        }
        names(ah.query)[ah.query.which[1]]
        } else {
        override.code
        }
    }
    annotation_key <- find.annotationhub.name(annotation_genus_species, annotation_key_override)
    orgdb <- ah[[annotation_key]]
    return(orgdb)
}

#' Plot an interactive PCA plot
#'
#' @param rld DESeqTransform object output by varianceStabilizingTransformation() or rlog()
#' @param intgroup character vector of names in colData(x) to use for grouping
#'
#' @return Handle to ggplot with added label field in aes_string() for plotting with ggplotly()
plotPCA.ly <- function(rld, intgroup){
  mat <- plotPCA(rld, intgroup, returnData=TRUE)
  pv <- attr(mat, 'percentVar')
  p <- ggplot(data=mat, aes_string(x='PC1', y='PC2', color='group', label='name')) +
          geom_point(size=3) + xlab(paste0('PC1: ', round(pv[1]*100), '% variance')) +
          ylab(paste0('PC2: ', round(pv[2]*100), '% variance')) + coord_fixed()
  return(p)
}

#' Plot a MA plot labeled with selected genes
#'
#' @param res.list data.frame (or list) with log2FoldChange and padj columns (or elements).
#'        Row names of the data.frame should be gene IDs/symbols and match the format of lab.genes
#' @param fdr.thres FDR threshold for defining statistical significance
#' @param fc.thres log2FoldChange cutoff for defining statistical significance
#' @param fc.lim User-defined limits for y-axis (log2FC). If NULL, this is defined as
#'        the (floor, ceiling) of range(res$log2FoldChange)
#' @param genes.to.label Genes to label on the MA plot. NULL, by default
#' @param col Column of `res` in which to look for `genes.to.label`. If NULL,
#'        rownames are used.
#'
#' @return Handle to ggplot
plotMA.label <- function(res,
                         fdr.thres=0.1,
                         fc.thres=0,
                         fc.lim=NULL,
                         genes.to.label=NULL,
                         col=NULL
                         ){
  # TODO: Add ggrastr option for points
  genes.to.label <- as.character(genes.to.label)
  nna <- sum(is.na(genes.to.label))
  if (nna > 0){
      warning(paste("Removing", nna, "NAs from gene list"))
      genes.to.label <- genes.to.label[!is.na(genes.to.label)]
  }
  # convert res to data frame
  res <- data.frame(res)

  # if y limits not specified
  if(is.null(fc.lim)){
    fc.lim <- range(res$log2FoldChange, na.rm=TRUE)
    fc.lim[1] <- floor(fc.lim[1])
    fc.lim[2] <- ceiling(fc.lim[2])
  }

  # get data frame of genes outside plot limits
  up.max <- res[res$log2FoldChange > fc.lim[2],]
  up.max$log2FoldChange <- rep(fc.lim[2], dim(up.max)[1])
  up.max <- data.frame(genes=rownames(up.max), up.max)

  down.max <- res[res$log2FoldChange < fc.lim[1],]
  down.max$log2FoldChange <- rep(fc.lim[1], dim(down.max)[1])
  down.max <- data.frame(genes=rownames(down.max), down.max)

  # get data frame of DE genes
  de.list <- res[res$padj < fdr.thres &
                 !is.na(res$padj) &
                 abs(res$log2FoldChange) >= fc.thres,]
  de.list <- data.frame(genes=rownames(de.list), de.list)

  # get data frame of DE genes outside plot limits
  up.max.de <- up.max[rownames(up.max) %in% rownames(de.list),]
  down.max.de <- down.max[rownames(down.max) %in% rownames(de.list),]

  # create ggplot with appropriate layers
  p <- ggplot(res, aes(baseMean, log2FoldChange)) +
    geom_point(col="gray40") + scale_x_log10() + ylim(fc.lim[1], fc.lim[2]) +
    theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

  p <- p + geom_hline(yintercept = 0, col="red", size=2, alpha=0.5)     # add horizontal line
  p <- p + geom_point(data=up.max, col="gray40", pch=2)                 # add points above max y
  p <- p + geom_point(data=down.max, col="gray40", pch=6)               # add points below min y

  p <- p + geom_point(data=de.list, col="red")                          # add DE points
  p <- p + geom_point(data=up.max.de, col="red", pch=2)                 # add DE points above max y
  p <- p + geom_point(data=down.max.de, col="red", pch=6)               # add DE points below min y


  if(!is.null(genes.to.label)){
    # get data frame of genes to be labeled
    if (!is.null(col)){
        res$gene.labels <- res[,col]
    } else {
        res$gene.labels <- rownames(res)
    }

    label.list <- res[res$gene.labels %in% genes.to.label,]
    #label.list <- data.frame(genes=rownames(label.list), label.list)

    # label genes outside limits
    up.max.idx <- rownames(label.list) %in% rownames(up.max)
    down.max.idx <- rownames(label.list) %in% rownames(down.max)

    if(sum(up.max.idx) > 0){
      label.list$log2FoldChange[up.max.idx] <- rep(fc.lim[2], sum(up.max.idx))
    }

    if(sum(down.max.idx) > 0){
      label.list$log2FoldChange[down.max.idx] <- rep(fc.lim[1], sum(down.max.idx))
    }

    # add labels
    p <- p + geom_point(data=label.list, col="black", pch=1, size=3)
    p <- p + geom_label_repel(data=label.list, aes(label=label.list$gene.labels, fontface="italic"))
  }
  return(p)

}

#' Plot a clustered heatmap of samples
#'
#' @param rld DESeqTransform object, typically output from running rlog()
#' @param colData Dataframe of metadata, used for annotating heatmap
#' @param cols.for.grouping Columns in colData to annotate heatmap with
#'
#' @return matrix of sample distances
plot.heatmap <- function(rld, colData, cols.for.grouping){
    sampleDists <- dist(t(assay(rld)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- colnames(rld)
    colors <- colorRampPalette(rev(brewer.pal(9, 'Blues')))(255)
    df <- as.data.frame(colData(rld)[, cols.for.grouping])
    colnames(df) <- cols.for.grouping
    rownames(df) <- colnames(rld)
    heatmaply(sampleDistMatrix,
              scale='none',
              col=colors,
              row_side_colors=df,
              showticklabels=c(FALSE,TRUE))
}


#' Plot heatmap of most varying genes
#'
#' @param rld DESeqTransform object, typically output from running rlog()
#' @param colData Dataframe of metadata, used for annotating heatmap
#' @param n Number of genes to include
#'
#' @return Side effect is to plot the heatmap
vargenes.heatmap <- function(rld, cols.for.grouping, n=50){
  # TODO: heatmaply?
  # TODO: allow alternative gene IDs
  topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE), n)
  mat <- assay(rld)[topVarGenes,]
  mat <- mat - rowMeans(mat)
  df <- as.data.frame(colData(rld)[, cols.for.grouping])
  rownames(df) <- colnames(rld)
  colnames(df) <- cols.for.grouping
  pheatmap(mat, annotation_col=df, cluster_cols=TRUE)
}


#' Simple wrapper of cat() that makes markdown text easier to print
#'
#' @param ... Arguments to print
#'
#' Make sure you're in an R code chunk with `results="asis"` set.
mdcat <- function(...){
  cat('\n\n', ..., ' \n\n', sep='', fill=1500)
}


#' Convert a vector of BAM filenames into a vector of samplenames.
#'
#' @param x Character vector of BAM filenames (e.g., from the header of
#'        featureCounts)
#'
lcdbwf.samplename <- function(x) {
    x <- x %>%
        str_remove_all('data/rnaseq_samples/') %>%
        str_remove_all('.cutadapt.bam') %>%
        str_split(fixed('/'), simplify=TRUE)
    x[,1]
}


#' Load a single combined featureCounts table into a DESeq object.
#'
#' @param filename Filename containing featureCounts combined output
#' @param sampletable data.frame containing at least sample names as first
#'        column
#' @param design Model used for creating DESeq object
#' @param sample.func Function that will be applied to each column name
#'        to align it to the input sampletable. Takes a character vector as
#'        input and returns a character vector of adjusted strings.
#' @param subset.counts If TRUE, then counts are subsetted by the provided
#'        sampletable. That is, only the counts columns that match a rowname of
#'        the provided sampletable are included. If FALSE (default), an error
#'        is raised reporting the differences in sampletable and counts data.
#'        It is always an error if counts data does not have an entry for
#'        a sample in the sampletable.
#'
#' @return DESeq object
#'
#' Additional args are passed to DESeq2::DESeqDataSetFromMatrix.
DESeqDataSetFromCombinedFeatureCounts <- function(filename, sampletable, design, sample.func=lcdbwf.samplename, subset.counts=FALSE, ...){

    # The sampletable may be data.frame or tibble; if it's a tibble then it
    # likely doesn't have rownames. So in this function we assume that it's the
    # first column that contains the samplenames.

    # Read in the counts TSV, use gene ID as rownames, and get rid of the Chr,
    # Start, End, Strand, Length columns
    m <- read_tsv(filename, comment="#") %>%
        remove_rownames %>%
        column_to_rownames('Geneid') %>%
        dplyr::select(-(1:5)) %>%
        as.data.frame

    # The column names of the imported table are not guaranteed to match the
    # samples in the metadata, so we need to make sure things match up
    # correctly when renaming columns.
    #
    # sample.func should convert filenames (which are columns in featurecounts
    # table) with samples (as listed in the sampletable)
    x <- colnames(m) %>% sample.func

    samplenames <- sampletable[,1]
    counts.not.sampletable <- setdiff(x, samplenames)
    sampletable.not.counts <- setdiff(samplenames, x)


    if (!all(x %in% samplenames)){
        # sampletable is missing:
        if (!subset.counts){
            stop(paste('The following samples are in the counts data but not the sampletable. If this is intended, consider using `subset.counts=TRUE` to remove them from the counts:', paste(counts.not.sampletable, collapse=', ')))
        }
    }
    if (!all(samplenames %in% x)){
        stop(
           paste(
             'The following samples are in the sampletable but not in the counts data. Check sample.func?',
             paste(sampletable.not.counts, collapse=', '))
        )
    }

    colnames(m) <- x

    # We should be good -- so reorder columns in `m` to match sampletable
    m <- m[, sampletable[,1]]

    # Before creating the dds object, we will drop levels from factors, in case
    # that was not done before sending in the (possibly subsetted) sampletable.

    sampletable <- droplevels(as.data.frame(sampletable))

    object <- DESeqDataSetFromMatrix(countData=m, colData=sampletable, design=design, ...)
    return(object)
}

#' Load featureCounts output into a DESeq object.
#'
#' Revised version of DESeq2::DESeqDataSetFromHTSeqCount to handle the
#' featureCounts default output format, which contains many more columns.
#'
#' @param sampleTable data.frame containing at least "featurecounts.path" column
#' @param directory Paths to featureCounts output are relative to this dir
#' @param design Model used for creating DESeq object
#'
#' @return DESeq object
#'
#' Additional args are passed to DESeq2::DESeqDataSetFromMatrix.
DESeqDataSetFromFeatureCounts <- function (sampleTable, directory='.', design,
                                           ignoreRank=FALSE,  ...)
{
  l <- lapply(
    as.character(sampleTable[, 'featurecounts.path']),
    function(fn) read.table(file.path(directory, fn), stringsAsFactors=FALSE, skip=2)
  )
  if (!all(sapply(l, function(a) all(a$V1 == l[[1]]$V1))))
    stop("Gene IDs in first column differ between files")
  tbl <- sapply(l, function(a) a$V7)
  colnames(tbl) <- sampleTable[, 1]
  rownames(tbl) <- l[[1]]$V1
  rownames(sampleTable) <- sampleTable[, 1]
  object <- DESeqDataSetFromMatrix(countData=tbl, colData=sampleTable[, -grepl('path', colnames(sampleTable)),
                                   drop=FALSE], design=design, ignoreRank, ...)
  return(object)
}

#' Load Salmon quantification data into a DESeq object
#'
#' @param sampleTable data.frame containing at least "salmon.path" column
#' @param design Model used for creating DESeq object
#'
#' @return DESeq object with transcript-level counts
#'
#' Additional args are passed to DESeq2::DESeqDataSetFromMatrix.
DESeqDataSetFromSalmon <- function (sampleTable, design,
                                           ignoreRank=FALSE,  ...)
{
    txi <- tximport(sampleTable[, 'salmon.path'], type='salmon', txOut=TRUE)
    object <- DESeqDataSetFromTximport(txi, colData=sampleTable[, -grepl('path', colnames(sampleTable)),
                                       drop=FALSE], design=design, ignoreRank, ...)
    return(object)
}

#' Make a single dds object
#'
#' This function is intended to be applied to a list of design data.
#'
#' @param design_data Named list of 2 to 4 items. At least "sampletable" (the
#'           colData, as a data.frame or tibble) and "design" (e.g., `~group`)
#'           are required. The optional named items are "file", which is the
#'           featureCounts file containing counts for all samples, and "args"
#'           which is a list of arguments to be passed to the constructor
#'           (e.g., `args=list(subset.counts=TRUE))`.
#' @param salmon.files If you want the dds objects to be generated from salmon
#'           counts, then provide the txi object from tximport. Otherwise,
#'           leave as NULL to use featureCounts. The value of this argument is
#'           used for all dds objects in the returned list (that is, the
#'           returned list cannot have a mix of salmon and featureCounts; if
#'           you want both you'll need to call this function twice).
#' @param combine.by The column to collapse technical replicates by. Rows in
#'           the sampletable that share the same value of this column will be
#'           combined using DESeq2::collapseReplicates.
#' @param ... Additional arguments will be passed on to the DESeq() call (e.g.,
#'           parallel, fitType, etc)
#' @param remove.version If TRUE, gene (or transcript) version information --
#'           the ".1" in "ENSG0000102345.1" -- will be stripped off.
make.dds <- function(design_data, salmon.files=NULL, combine.by=NULL,
                     remove.version=FALSE, ...){
    colData <- pluck(design_data, 'sampletable')
    design <- pluck(design_data, 'design')
    location <- pluck(design_data, 'file',
                      .default='../data/rnaseq_aggregation/featurecounts.txt')
    arg_list <- pluck(design_data, 'args')

    if (is.null(salmon.files)) {
        dds <- exec(
            DESeqDataSetFromCombinedFeatureCounts,
                location,
                sampletable=colData,
                design=design,
                !!!arg_list)
    } else {
        dds <- exec(
            DESeqDataSetFromTximport,
            salmon.files,
            colData=colData[, -grepl('path', colnames(colData)), drop=FALSE],
            design=design,
            !!!arg_list)
    }

    if (remove.version){
        rownames(dds) <- sapply(strsplit(rownames(dds), '.', fixed=TRUE),
                                function (x) {ifelse(grepl('_', x[2]),
                                              paste(x[1], x[2], sep='.'),
                                              x[1])}
                                )
    }

    if(!is.null(combine.by)){
        dds <-collapseReplicates(dds, dds[[combine.by]])
    }

    dds <- DESeq(dds, ...)
    return(dds)
}

#' Make a list of dds objects
#'
#' Helper function to construct a list of dds objects. The `make.dds` function
#' does all the work; this just sets up some sane defaults and does the map()
#' call.
#'
#' @param deseq_obj_list A named list of lists. Each list is used as the first
#'           argument to `make.dds`; see the documentation of that function for
#'           details.
#'
#' @return A list of dds objects.
#'
make.dds.list <- function(deseq_obj_list, salmon.files=NULL, combine.by=FALSE,
                          remove.version=TRUE, ...){
    dds_list <- map(deseq_obj_list, make.dds, salmon.files, combine.by,
                    remove.version, ...)
    return(dds_list)
}

#' Compute label for one component of an arbitrary design matrix
#'
#' @param pattern.mat single row from unique(model.matrix) call from a dds object
#'
#' @return formatted column name

format.group.name <- function(pattern.mat) {
    my.res <- "counts"
    # enforce input format
    stopifnot(is.vector(pattern.mat))
    # if there's only one element
    if (length(pattern.mat) == 1) {
        # it must be the intercept
    my.res <- paste(my.res, "groupMean", sep=".")
    } else { # more than one element
        # for every non-intercept entry
        for (i in 2:length(pattern.mat)) {
        # if 1, include
        if (pattern.mat[i] == 1) my.res <- paste(my.res, names(pattern.mat[i]), sep=".")
        # otherwise don't
        else my.res <- paste(my.res, ".not", names(pattern.mat[i]), sep="")
    }
    }
    my.res
}

#' Plot a gene's normalized counts across samples
#'
#' @param gene Gene ID
#' @param dds DESeq object from which to extract counts
#'
#' @return ggplot object
my.counts <- function(gene, dds, label=NULL, intgroup='group'){

  # Assumption: color genes by group
  geneCounts <- plotCounts(dds, gene=gene, intgroup=intgroup, returnData=TRUE)
  p <- ggplot(geneCounts, aes_string(x=intgroup, y='count', color=intgroup, group=intgroup)) +
    scale_y_log10() +
    geom_point(position=position_jitter(width=.1, height=0),  size=3) +
    geom_line(color='#000000') +
    ggtitle(gene)

  if (!is.null(label)){
    p <- p + ggtitle(label)
  }
  return(p)
}

#' Re-order results by log2FoldChange
#'
#' @param res DESeq2 results object
#' @param reverse If TRUE then sort high-to-low
#'
#' @return Re-ordered results object
lfc.order <- function(res, reverse=FALSE){
    res.na <- res[!is.na(res$log2FoldChange),]
    if (!reverse){
        return(res.na[order(res.na$log2FoldChange),])
    }
    if (reverse){
        return(res.na[rev(order(res.na$log2FoldChange)),])
    }
}

#' Re-order results by adjusted pval
#'
#' @param res DESeq2 results object
#' @param reverse If TRUE then sort high-to-low
#'
#' @return Re-ordered results object
padj.order <- function(res, reverse=FALSE){
    res.na <- res[!is.na(res$log2FoldChange) & !is.na(res$log2FoldChange),]
    if (!reverse){
        res.na <- res.na[res.na$log2FoldChange > 0,]
    } else {
        res.na <- res.na[res.na$log2FoldChange < 0,]
    }
    return(res.na[order(res.na$padj),])
}


#' Plot normalized gene counts for the top N genes
#'
#' @param sorted DESeq2 results object
#' @param n Number of top genes to plot
#' @param func Plotting function to call on each gene. Signature is func(gene,
#'             dds, label=label, ...)
#' @param dds DESeq2 object
#' @param label Vector of column names in `res` from which to add a label to
#'        the gene (e.g., c('symbol', 'alias'))
#' @param ... Additional arguments that are passed to `func`
#' @return Side effect is to create plot
top.plots <- function(res, n, func, dds, add_cols=NULL, ...){
    ps <- list()
    for (i in seq(n)){
        gene <- rownames(res)[i]
        add_label <- as.character(as.data.frame(res)[i, add_cols])
        add_label <- add_label[!is.na(add_label)]
        label <- paste(gene, add_label, sep=' | ')
        if (length(label) == 0){
          label <- NULL
        }
        ps[[gene]] <- func(gene, dds, label=label, ...)
    }
    grid.arrange(grobs=ps)
}

#' Make dataframe of normalized gene counts
#'
#' @param dds DESeq2 object
#' @param res DESeq2 results object
#' @param sel.genes list of genes to consider
#' @param label column(s) to be included to the plot labels
#' @param pc count number to be added to the normalized counts
#'        Typically a pc of 0.5 is added to allow plotting in log scale
#'
#' @return dataframe
counts.df <- function(dds, res, sel.genes=NULL, label=NULL, rank.col='padj', pc=0.5) {
    # getting normalized counts copied from plotCounts()
    cnts <- counts(dds, normalized = TRUE, replaced = FALSE)
    # add 0.5 like plotCounts to plot 0 on log scale
    cnts <- cnts + pc
    # merge with res.i and colData
    df <- as.data.frame(cbind(cnts, res))
    # add label for plotting
    df <- df %>% mutate(label = paste(gene, !!!syms(label), sep=' | '))
    # subset to sel.genes if not NULL (then keep all)
    if (!is.null(sel.genes)) {
        df <- df %>%
            filter(gene %in% sel.genes)
    }
    # add rank
    df <- df %>%
        mutate(rank = rank(!!sym(rank.col), ties.method='first', na.last='keep')) %>%
        pivot_longer(colnames(cnts), names_to='samplename', values_to='normalized_counts')
    # add colData
    df <- merge(df, colData(dds), by='samplename') %>%
        as.data.frame()
    return(df)
}

#' Filter results by log2FoldChange sign
#'
#' @param res DESeq2 results object
#' @param reverse If TRUE then filter negative lfc
#'
#' @return filtered genes
lfc.filter <- function(res, reverse=FALSE){
    res.na <- res[!is.na(res$log2FoldChange),]
    if (!reverse){
        res.na <- res.na[res.na$log2FoldChange > 0,]
    } else {
        res.na <- res.na[res.na$log2FoldChange < 0,]
    }
    return(res.na[,'gene'])
}

#' Plot genes' normalized counts across samples
#'
#' @param df dataframe from which to extract genes normalized_counts
#' @param rank.nb rank number less or equal used to filter the genes to be plotted
#' @param no.aes whether to include the default aes or return the ggplot object
#' without aes
#' @param facet name of column to use for faceting
#'
#' @return ggplot object
counts.plot <- function(df, rank.nb=NULL, no.aes=FALSE, facet='label') {
    if (!is.null(rank.nb)) {
        df <- df %>%
            filter(rank <= rank.nb)
    }
    df <- df %>%
        arrange(rank) %>%
        mutate(facet = factor(!!!syms(facet), levels = unique(!!!syms(facet))))
    plt <- ggplot(df) +
        scale_y_log10() +
        geom_point(position=position_jitter(width=.1, height=0),  size=3) +
        geom_line(color='#000000') +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        facet_wrap(.~facet, ncol=1, scales='free_y')
    if (!no.aes) {
        plt <- plt +
            aes(y=normalized_counts, x=group, color=group)
    }
    return(plt)
}

#' Plot a histogram of raw pvals
#'
#' @param res DESeq2 results object
#'
#' @return Side effect is to create plot
pval.hist <- function(res){
    hist(res$pvalue[res$baseMean>1], breaks=0:20/20, col='grey50',
         border='white', xlab='P-value', main='Distribution of p-values')
}

#' Summarize DESeq2 results into a dataframe
#'
#' summary(res) prints out info; this function captures it into a dataframe
#'
#' @param res DESeq2 results object
#' @param dds DEseq2 object
#' @param alpha Alpha level at which to call significantly changing genes
#'
#' @return Dataframe of summarized results
my.summary <- function(res, dds, alpha, lfc.thresh=0, ...){
   if (missing(alpha)){
       alpha <- if (is.null(metadata(res)$alpha)){ 0.1 } else { metadata(res)$alpha }
   }
   notallzero <- sum(res$baseMean > 0)
   up <- sum(res$padj < alpha & res$log2FoldChange > lfc.thresh, na.rm=TRUE)
   down <- sum(res$padj < alpha & res$log2FoldChange < -lfc.thresh, na.rm=TRUE)
   filt <- sum(!is.na(res$pvalue) & is.na(res$padj))
   outlier <- sum(res$baseMean > 0 & is.na(res$pvalue))
   ft <- if(is.null(metadata(res)$filterThreshold)){ 0 } else { round(metadata(res)$filterThreshold) }
   # adjust width.cutoff as newline insertion causes this to return a df with
   # multiple duplicate rows!
   df <- data.frame(
                    total.annotated.genes=nrow(res),
                    total.nonzero.read.count=notallzero,
                    alpha=alpha,
                    lfcThreshold=lfc.thresh,
                    up=up,
                    down=down,
                    outliers=outlier,
                    low.counts=filt,
                    design=deparse(design(dds), width.cutoff=500L)
                    )
   return(df)
}

#' Convert "1/100" to 0.01.
#'
#' clusterProfiler report columns that are strings of numbers; this converts to
#' a fraction
#'
#' @param x Character vector to convert
get.frac <- function(x){
    y <- as.numeric(strsplit(x, '/')[[1]])
    return(y[[1]] / y[[2]])
}

#' Summarize and aggregate multiple GO results
#'
#' Convert a list of GO analysis results (one per ontology) to a single
#' dataframe, with additional label column filled in with `label` and with
#' a fraction column.
#'
#' @param ego List of clusterProfiler results objects
#' @param labels List of labels. For each name, a new column will be added and
#"        its contents will be the value
#'
#' @return dataframe
summarize.go <- function(ego, labels){
  lst <- list()
  for (name in names(ego)){
    d <- as.data.frame(ego[[name]])
    if (nrow(d) > 0){
      d$ontology <- name
      for (label in names(labels)){
        d[label] <- labels[[label]]
      }
      d$frac <- sapply(d$GeneRatio, get.frac)
    }
    lst[[name]] <- d
  }
  df <- do.call(rbind, lst)
  return(df)
}


#' Summarize KEGG results
#'
#' Summarize KEGG results and add `frac` and `label` columns
#'
#' @param ekg Results from clusterProfiler KEGG enrichment
#' @param label Attach this label to the "label" column
#'
#' @return Dataframe
summarize.kegg <- function(ekg, labels){
  d <- as.data.frame(ekg)
  if (nrow(d) > 0){
    d$ontology <- 'kegg'
    for (label in names(labels)){
      d[label] <- labels[[label]]
    }
    d$frac <- sapply(d$GeneRatio, get.frac)
  }
  return(d)
}

#' Attach additional information from an OrgDb to a results dataframe
#'
#' @param res DESeq2 results object
#' @param keytype Key type for rownames(res)
#' @param columns Character vector of columns to add from OrgDb to results
attach.info <- function(res, keytype='ENSEMBL', columns=c('SYMBOL', 'UNIPROT', 'ALIAS')){
    keys <- rownames(res)
    for (column in columns){
      label <- tolower(column)
      res[label] <- mapIds(orgdb, keys=keys, column=column, keytype=keytype, multiVal='first')
      res[[label]] <- ifelse(is.na(res[[label]]), keys, res[[label]])
    }
    # Put "gene" column as the first
    cn <- colnames(res)
    res$gene <- rownames(res)
    res <- res[, c('gene', cn)]
    return(res)
}


#' Combine everything in the results list into a single table
#'
#' @param res.list Named list of lists, where each sublist contains the following
#'                 names: c('res', 'dds', 'label'). "res" is a DESeqResults object,
#'                 "dds" is either the indexing label for the dds.list object or
#'                  the DESeq object, and "label" is a nicer-looking
#'                 label to use. NOTE: backwards compatibility with older versions
#'                  of lcdb-wf depends on no dds.list object being passed.
#'
#' @return Dataframe
summarize.res.list <- function(res.list, alpha, lfc.thresh, dds.list=NULL){
    slist <- list()
    for (name in names(res.list)){
        if(!is.null(dds.list)){
            x <- my.summary(res.list[[name]][['res']], dds.list[[ res.list[[name]][['dds']] ]], alpha, lfc.thresh)
        } else { x <- my.summary(res.list[[name]][['res']], res.list[[name]][['dds']], alpha, lfc.thresh)
        }
        rownames(x) <- res.list[[name]][['label']]
        slist[[name]] <- x
    }
    slist <- do.call(rbind, slist)
    return(slist)
}


#' Return index of up/down/changed genes
#'
#' @param x DESeq2 results object, or data.frame created from one
#' @param direction Direction in 'up', 'dn', 'down', 'ch', 'changed'
#' @param alpha FDR lower than this will be considered significant
#' @param thresh Log2 fold change threshold. If e.g. 2, will return < -2 and/or > 2, depending on the value of "direction"
#' @param return.names If TRUE, returns the rownames of selected genes; if FALSE return boolean index of length(x)
#'
#' @return Character vector of rownames (if return.names=TRUE) or boolean vector of genes selected.
get.sig <- function(x, direction='up', alpha=0.1, lfc.thresh=0, return.names=TRUE){
    if (direction == 'up'){
        idx <- (x$padj < alpha) & (x$log2FoldChange > lfc.thresh) & (!is.na(x$padj))
    } else if (direction %in% c('down', 'dn')){
        idx <- (x$padj < alpha) & (x$log2FoldChange < -lfc.thresh) & (!is.na(x$padj))
    } else if (direction %in% c('changed', 'ch')){
        idx <- (x$padj < alpha) & (abs(x$log2FoldChange) > lfc.thresh) & (!is.na(x$padj))
    }
    if (return.names){
        return(rownames(x)[idx])
    } else {
        return(idx)
    }
}

#' Prepare list for UpSet plots, but include rownames
#'
#' @param x List of sets to compare (same input as to UpSetR::fromList)
#'
#' @return data.frame of 1 and 0 showing which genes are in which sets
fromList.with.names <- function(lst){
    elements <- unique(unlist(lst))
    data <- unlist(lapply(lst, function(x) {
        x <- as.vector(match(elements, x))
    }))
    data[is.na(data)] <- as.integer(0)
    data[data != 0] <- as.integer(1)
    data <- data.frame(matrix(data, ncol = length(lst), byrow = F))

    # This is the only way in which this function differs from UpSetR::fromList
    rownames(data) <- elements

    data <- data[which(rowSums(data) != 0), ]
    names(data) <- names(lst)
    return(data)
}


#' Move rownames of data.frame to first column
#'
#' @param x data.frame
#' @param name Name to use as the new first column
#'
#' @return data.frame with original rownames moved to new first column with
#' optional name. rownames on the new data.frame are set to NULL.
rownames.first.col <- function(x, name='names'){
    orig.names <- colnames(x)
    x <- data.frame(rownames(x), x, row.names=NULL)
    colnames(x) <- c(name, orig.names)
    return(x)
}


#' lapply over a nested list
#'
#' @param x Nested list. Expected 2 levels deep, and each item is the same type
#' @param subfunc Function to apply over nested "leaf" items
#' @param ... Additional arguments to be passed to subfunc
#'
#' @return Nested list with same structure but tranformed values
nested.lapply <- function(x, subfunc, ...){
    lapply(x, lapply, subfunc, ...)
}

#' Split clusterProfiler output into one line per gene
#'
#' @param x Results from clusterProfiler. It is expected that the
#' clusterProfiler enrichment function was called with "readable=TRUE"
#'
#' @return data.frame with genes one per line instead of "/" separated in one
#' line. The rest of the original line is repeated for each gene.
#'
split.clusterProfiler.results <- function(x){
    df <- x@result
    # loop over all rows
    df.parse <- NULL
    for(k in 1:dim(df)[1]){
        g <- strsplit(as.character(df$geneID[k]), "/")[[1]]
        gg <- df[rep(k, length(g)),]
        gg$geneParse <- g
        if(is.null(df.parse)){
            df.parse <- gg
        } else {
            df.parse <- rbind(df.parse, gg)
        }
    }
    return(df.parse)
}


#' Writes out original and split clusterprofiler results
#'
#' @param res clusterProfiler results
#' @param cprof.folder Directory in which to save files. Will be created if needed.
#' @param Label to use for the results. Will generate a filename
#' cprof.folder/label.txt and cprof.folder/label_split.txt
#'
#' @return List of the two files that are created on disk.
write.clusterprofiler.results <- function(res, cprof.folder, label){
    dir.create(cprof.folder, showWarnings=FALSE, recursive=TRUE)
    filename.orig <- file.path(cprof.folder, paste0(label, '.txt'))
    write.table(res, file=filename.orig, sep='\t', quote=FALSE, row.names=FALSE)
    filename.split <- file.path(cprof.folder, paste0(label, '_split.txt'))
    res.split <- split.clusterProfiler.results(res)
    write.table(res.split, file=filename.split, sep='\t', quote=FALSE, row.names=FALSE)
    return(list(orig=filename.orig, split=filename.split))
}

