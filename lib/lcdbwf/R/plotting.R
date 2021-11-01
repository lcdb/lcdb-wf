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
#' @param genes_to_label genes_to_label on the MA plot. NULL, by default
#' @param col Column of `res` in which to look for `genes_to_label`. If NULL,
#'        rownames are used.
#'
#' @return Handle to ggplot
plotMA.label <- function(res,
                         fdr.thres=0.1,
                         fc.thres=0,
                         fc.lim=NULL,
                         genes_to_label=NULL,
                         label_column=NULL
                         ){
  # TODO: Add ggrastr option for points
  genes_to_label <- as.character(genes_to_label)
  nna <- sum(is.na(genes_to_label))
  if (nna > 0){
      warning(paste("Removing", nna, "NAs from gene list"))
      genes_to_label <- genes_to_label[!is.na(genes_to_label)]
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


  if(!is.null(genes_to_label)){
    # get data frame of genes to be labeled
    if (!is.null(label_column)){
        res$gene.labels <- res[,label_column]
    } else {
        res$gene.labels <- rownames(res)
    }

    label.list <- res[res$gene.labels %in% genes_to_label,]
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

#' Plot a volcano plot labeled with selected genes
#'
#' @param res.list data.frame (or list) with log2FoldChange and padj columns (or elements).
#'        Row names of the data.frame should be gene IDs/symbols and match the format of lab.genes
#' @param fdr.thres FDR threshold for defining statistical significance
#' @param fc.thres log2FoldChange cutoff for defining statistical significance
#' @param fc.lim User-defined limits for x-axis (log2FC). If NULL, this is defined as
#'        the (floor, ceiling) of range(res$log2FoldChange)
#' @param genes_to_label genes_to_label on the volcano plot. NULL, by default
#' @param col Column of `res` in which to look for `genes_to_label`. If NULL,
#'        rownames are used.
#'
#' @return Handle to ggplot
plot.volcano.label <- function(res,
                         fdr.thres=0.1,
                         fc.thres=0,
                         fc.lim=NULL,
                         genes_to_label=NULL,
                         label_column=NULL
                         ){
  # TODO: Add ggrastr option for points
  genes_to_label <- as.character(genes_to_label)
  nna <- sum(is.na(genes_to_label))
  if (nna > 0){
      warning(paste("Removing", nna, "NAs from gene list"))
      genes_to_label <- genes_to_label[!is.na(genes_to_label)]
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
  p <- ggplot(res, aes(log2FoldChange, -log10(padj))) +
    geom_point(col="gray40") + xlim(fc.lim[1], fc.lim[2]) +
    theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

  p <- p + geom_point(data=up.max, col="gray40", pch=2)                 # add points above max y
  p <- p + geom_point(data=down.max, col="gray40", pch=6)               # add points below min y

  p <- p + geom_point(data=de.list, col="red")                          # add DE points
  p <- p + geom_point(data=up.max.de, col="red", pch=2)                 # add DE points above max y
  p <- p + geom_point(data=down.max.de, col="red", pch=6)               # add DE points below min y


  if(!is.null(genes_to_label)){
    # get data frame of genes to be labeled
    if (!is.null(label_column)){
        if (!(label_column %in% colnames(res))){
            stop(paste(label_column, "is not a column in the results object; columns are:", colnames(res)))
        }
        res$gene.labels <- res[,label_column]
    } else {
        res$gene.labels <- rownames(res)
    }

    label.list <- res[res$gene.labels %in% genes_to_label,]
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

#' Barplot of size factors by sample
#'
#' @param dds DESeqDataSet object
sizefactors_barplot <- function(dds){
    dds <- DESeq2::estimateSizeFactors(dds)
    sf <- DESeq2::sizeFactors(dds)
    sf <- sf[order(sf)] %>%
            enframe(value = 'Size Factor')
    p <- ggplot(sf) +
        aes(x=reorder(name, `Size Factor`), y=`Size Factor`) +
        xlab('sample name') +
        geom_col() +
        theme_bw() +
        theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
    return(p)
}

#' Scatterplot of size factors vs total read count
#'
#' @param dds DESeqDataSet object
sizefactors_vs_total <- function(dds){
    dds <- DESeq2::estimateSizeFactors(dds)
    sf <- DESeq2::sizeFactors(dds)
    sf <- sf[order(sf)] %>%
            enframe(value = 'Size Factor')
    trc <- colSums(counts(dds)) %>%
            enframe(value = 'Total Read Count')
    trc_vs_sf <- full_join(sf, trc, by='name')
    p <- ggplot(data=trc_vs_sf, aes_string(x="`Total Read Count`", y="`Size Factor`", label='name')) +
        geom_point(size=3) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    return(p)
}
