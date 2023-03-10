#' Plot an interactive PCA plot
#'
#' @param rld DESeqTransform object output by varianceStabilizingTransformation() or rlog()
#' @param intgroup character vector of names in colData(x) to use for grouping
#'
#' @return Handle to ggplot with added label field in aes_string() for plotting with ggplotly()
plotPCA.ly <- function(rld, intgroup){
  mat <- DESeq2::plotPCA(rld, intgroup, returnData=TRUE)
  pv <- attr(mat, 'percentVar')
  p <- ggplot2::ggplot(data=mat, aes_string(x='PC1', y='PC2', color='group', label='name')) +
          geom_point(size=3) + xlab(paste0('PC1: ', round(pv[1]*100), '% variance')) +
          ylab(paste0('PC2: ', round(pv[2]*100), '% variance')) + coord_fixed()
  return(p)
}

#' Plot a MA plot labeled with selected genes
#'
#' @param res.list data.frame (or list) with log2FoldChange and padj columns (or elements).
#'        Row names of the data.frame should be gene IDs/symbols and match the format of lab.genes
#' @param fdr_thres fdr_threshold for defining statistical significance
#' @param fd_thres log2FoldChange cutoff for defining statistical significance
#' @param fc_lim User-defined limits for y-axis (log2FC). If NULL, this is defined as
#'        the (floor, ceiling) of range(res$log2FoldChange)
#' @param genes_to_label genes_to_label on the MA plot. NULL, by default
#' @param col Column of `res` in which to look for `genes_to_label`. If NULL,
#'        rownames are used.
#'
#' @return Handle to ggplot
plotMA_label <- function(res,
                         fdr_thres=0.1,
                         fd_thres=0,
                         fc_lim=NULL,
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
  if(is.null(fc_lim)){
    fc_lim <- range(res$log2FoldChange, na.rm=TRUE)
    fc_lim[1] <- floor(fc_lim[1])
    fc_lim[2] <- ceiling(fc_lim[2])
  }

  # get data frame of genes outside plot limits
  up.max <- res[res$log2FoldChange > fc_lim[2],]
  up.max$log2FoldChange <- rep(fc_lim[2], dim(up.max)[1])
  up.max <- data.frame(genes=rownames(up.max), up.max)

  down.max <- res[res$log2FoldChange < fc_lim[1],]
  down.max$log2FoldChange <- rep(fc_lim[1], dim(down.max)[1])
  down.max <- data.frame(genes=rownames(down.max), down.max)

  # get data frame of DE genes
  de.list <- res[res$padj < fdr_thres &
                 !is.na(res$padj) &
                 abs(res$log2FoldChange) >= fd_thres,]
  de.list <- data.frame(genes=rownames(de.list), de.list)

  # get data frame of DE genes outside plot limits
  up.max.de <- up.max[rownames(up.max) %in% rownames(de.list),]
  down.max.de <- down.max[rownames(down.max) %in% rownames(de.list),]

  # create ggplot with appropriate layers
  p <- ggplot2::ggplot(res, aes(baseMean, log2FoldChange)) +
    ggplot2::geom_point(col="gray40") + scale_x_log10() + ylim(fc_lim[1], fc_lim[2]) +
    theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

  p <- p + ggplot2::geom_hline(yintercept = 0, col="red", size=2, alpha=0.5)     # add horizontal line
  p <- p + ggplot2::geom_point(data=up.max, col="gray40", pch=2)                 # add points above max y
  p <- p + ggplot2::geom_point(data=down.max, col="gray40", pch=6)               # add points below min y

  p <- p + ggplot2::geom_point(data=de.list, col="red")                          # add DE points
  p <- p + ggplot2::geom_point(data=up.max.de, col="red", pch=2)                 # add DE points above max y
  p <- p + ggplot2::geom_point(data=down.max.de, col="red", pch=6)               # add DE points below min y


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
      label.list$log2FoldChange[up.max.idx] <- rep(fc_lim[2], sum(up.max.idx))
    }

    if(sum(down.max.idx) > 0){
      label.list$log2FoldChange[down.max.idx] <- rep(fc_lim[1], sum(down.max.idx))
    }

    # add labels
    p <- p + ggplot2::geom_point(data=label.list, col="black", pch=1, size=3)
    p <- p + ggrepel::geom_label_repel(data=label.list, aes(label=label.list$gene.labels, fontface="italic"))
  }
  return(p)

}

#' Plot a volcano plot labeled with selected genes
#'
#' @param res.list data.frame (or list) with log2FoldChange and padj columns (or elements).
#'        Row names of the data.frame should be gene IDs/symbols and match the format of lab.genes
#' @param fdr_thres fdr_threshold for defining statistical significance
#' @param fd_thres log2FoldChange cutoff for defining statistical significance
#' @param fc_lim User-defined limits for x-axis (log2FC). If NULL, this is defined as
#'        the (floor, ceiling) of range(res$log2FoldChange)
#' @param genes_to_label genes_to_label on the volcano plot. NULL, by default
#' @param col Column of `res` in which to look for `genes_to_label`. If NULL,
#'        rownames are used.
#'
#' @return Handle to ggplot
plot_volcano_label <- function(res,
                         fdr_thres=0.1,
                         fd_thres=0,
                         fc_lim=NULL,
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
  if(is.null(fc_lim)){
    fc_lim <- range(res$log2FoldChange, na.rm=TRUE)
    fc_lim[1] <- floor(fc_lim[1])
    fc_lim[2] <- ceiling(fc_lim[2])
  }

  # get data frame of genes outside plot limits
  up.max <- res[res$log2FoldChange > fc_lim[2],]
  up.max$log2FoldChange <- rep(fc_lim[2], dim(up.max)[1])
  up.max <- data.frame(genes=rownames(up.max), up.max)

  down.max <- res[res$log2FoldChange < fc_lim[1],]
  down.max$log2FoldChange <- rep(fc_lim[1], dim(down.max)[1])
  down.max <- data.frame(genes=rownames(down.max), down.max)

  # get data frame of DE genes
  de.list <- res[res$padj < fdr_thres &
                 !is.na(res$padj) &
                 abs(res$log2FoldChange) >= fd_thres,]
  de.list <- data.frame(genes=rownames(de.list), de.list)

  # get data frame of DE genes outside plot limits
  up.max.de <- up.max[rownames(up.max) %in% rownames(de.list),]
  down.max.de <- down.max[rownames(down.max) %in% rownames(de.list),]

  # create ggplot with appropriate layers
  p <- ggplot(res, aes(log2FoldChange, -log10(padj))) +
    geom_point(col="gray40") + xlim(fc_lim[1], fc_lim[2]) +
    theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

  p <- p + ggplot2::geom_point(data=up.max, col="gray40", pch=2)                 # add points above max y
  p <- p + ggplot2::geom_point(data=down.max, col="gray40", pch=6)               # add points below min y

  p <- p + ggplot2::geom_point(data=de.list, col="red")                          # add DE points
  p <- p + ggplot2::geom_point(data=up.max.de, col="red", pch=2)                 # add DE points above max y
  p <- p + ggplot2::geom_point(data=down.max.de, col="red", pch=6)               # add DE points below min y


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
      label.list$log2FoldChange[up.max.idx] <- rep(fc_lim[2], sum(up.max.idx))
    }

    if(sum(down.max.idx) > 0){
      label.list$log2FoldChange[down.max.idx] <- rep(fc_lim[1], sum(down.max.idx))
    }

    # add labels
    p <- p + ggplot2::geom_point(data=label.list, col="black", pch=1, size=3)
    p <- p + ggrepel::geom_label_repel(data=label.list, aes(label=label.list$gene.labels, fontface="italic"))
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
plot_heatmap <- function(rld, colData, cols_for_grouping){
    sampleDists <- dist(t(assay(rld)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- colnames(rld)
    colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'Blues')))(255)
    df <- as.data.frame(colData(rld)[, cols_for_grouping])
    colnames(df) <- cols_for_grouping
    rownames(df) <- colnames(rld)
    heatmaply::heatmaply(sampleDistMatrix,
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
vargenes_heatmap <- function(rld, cols_for_grouping, n=50){
  # TODO: heatmaply?
  # TODO: allow alternative gene IDs
  topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE), n)
  mat <- assay(rld)[topVarGenes,]
  mat <- mat - rowMeans(mat)
  df <- as.data.frame(colData(rld)[, cols_for_grouping])
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
  p <- ggplot2::ggplot(geneCounts, aes_string(x=intgroup, y='count', color=intgroup, group=intgroup)) +
    scale_y_log10() +
    ggplot2::geom_point(position=position_jitter(width=.1, height=0),  size=3) +
    ggplot2::geom_line(color='#000000') +
    ggplot2::ggtitle(gene)

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
    gridExtra::grid.arrange(grobs=ps)
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
        dplyr::arrange(rank) %>%
        dplyr::mutate(facet = factor(!!!syms(facet), levels = unique(!!!syms(facet))))
    plt <- ggplot2::ggplot(df) +
        ggplot2::scale_y_log10() +
        ggplot2::geom_point(position=position_jitter(width=.1, height=0),  size=3) +
        ggplot2::geom_line(color='#000000') +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ggplot2::facet_wrap(.~facet, ncol=1, scales='free_y')
    if (!no.aes) {
        plt <- plt +
            ggplot2::aes(y=normalized_counts, x=group, color=group)
    }
    return(plt)
}

#' Plot a histogram of raw pvals
#'
#' This is edited from the DESeq2 vignette, from the section about independent
#' filtering. The resulting histogram indicates pvals for those genes kept and
#' removed before multiple testing adjustment.
#'
#' @param res DESeq2 results object
#'
#' @return Side effect is to create plot
pval_hist <- function(res){
  use <- res$baseMean > metadata(res)$filterThreshold
  h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
  h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
  colori <- c(`counts too low`='khaki', `pass`="powderblue")
  df <- rbind(data.frame(x=h1$mids, counts=h1$counts, label='counts too low'),
              data.frame(x=h2$mids, counts=h2$counts, label='pass')
              )
  plt <- ggplot2::ggplot(df, aes(x=x, y=counts, fill=label)) +
    geom_bar(stat = 'identity', color='gray20') +
    theme_classic() +
    scale_fill_manual(values=c("#EBE379", "#A3DAE0")) +
    xlab('p-value') +
    ylab('frequency') +
    theme(legend.position = c(0.8, 0.8))
  return(plt)
}

#' Barplot of size factors by sample
#'
#' @param dds DESeqDataSet object
sizefactors_barplot <- function(dds){
    dds <- DESeq2::estimateSizeFactors(dds)
    sf <- DESeq2::sizeFactors(dds)
    sf <- sf[order(sf)] %>%
            tibble::enframe(value = 'Size Factor')
    p <- ggplot2::ggplot(sf) +
        ggplot2::aes(x=reorder(name, `Size Factor`), y=`Size Factor`) +
        ggplot2::xlab('sample name') +
        ggplot2::geom_col() +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, vjust=0.5, hjust=1))
    return(p)
}

#' Scatterplot of size factors vs total read count
#'
#' @param dds DESeqDataSet object
sizefactors_vs_total <- function(dds){
    dds <- DESeq2::estimateSizeFactors(dds)
    sf <- DESeq2::sizeFactors(dds)
    sf <- sf[order(sf)] %>%
            tibble::enframe(value = 'Size Factor')
    trc <- colSums(counts(dds)) %>%
            tibble::enframe(value = 'Total Read Count')
    trc_vs_sf <- dplyr::full_join(sf, trc, by='name')
    p <- ggplot2::ggplot(data=trc_vs_sf, aes_string(x="`Total Read Count`", y="`Size Factor`", label='name')) +
        ggplot2::geom_point(size=3) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
    return(p)
}


#' Version of DESeq2::plotSparsity with ggplot2
plotSparsity2 <- function(dds){
  x <- counts(dds, normalize=TRUE)
  metadata <- S4Vectors::elementMetadata(dds) %>% as.data.frame()
  metadata$rowsums <- rowSums(x)
  metadata$rowmax <- apply(x, 1, max)

  ggplot2::ggplot(
    metadata %>%
      filter(rowsums > 0)) +
    ggplot2::aes(x=log10(rowsums), y=rowmax/rowsums) +
    ggplot2::geom_point() +
    ggplot2::xlab('sum of counts per gene') +
    ggplot2::ylab('max count / sum')
}




#' Plot a scatterplot of two contrasts' LFCs, color-coded by significance
#'
#' This is edited from the DESeq2 vignette, from the section about independent
#' filtering. The resulting histogram indicates pvals for those genes kept and
#' removed before multiple testing adjustment.
#'
#' @param res_i DESeq2 results object
#' @param res_j second DESeq2 results object
#' @param padj.thr float, p.adj threshold
#' @param name.col string, gene name column to merge the 2 results, also used for labelling plots
#' @param label_i, label_j string, label for res_i and res_j
#' @param return.values boolean,  whether to return the ggplot object (FALSE) or the dataframe (TRUE)
#' @param' color.palette list of string, colors to use for significance categories 'Both - same LFC sign',
#'         'Both - opposite LFC sign 'None', label_i, label_j

#' @return Either returns a ggplot object of the scatterplot, or the corresponding dataframe if return.values=TRUE

lfc_scatter <- function(res_i, res_j, padj.thr=0.1, name.col='SYMBOL', label_i=NULL, label_j=NULL,
                        return.values=FALSE, color.palette=c('#FF3333', "#FF6699", '#999999', '#66CCCC', '#0072B2')) {
    #  colors from color-blind palette                        red        pink       grey       cyan       blue
    # check whether the genes match in res_i and res_j, emits a warning if not
    diff.genes <- c(setdiff(rownames(res_i), rownames(res_j)),
                    setdiff(rownames(res_j), rownames(res_i))) %>%
                    unlist() %>%
                    unique()
    if( length(diff.genes) > 0 ) {
        warning(paste0(length(diff.genes),
                       ' genes were discarded because found in one res but not the other'))
    }

    # use generic labels if not provided
    if (is.null(label_i)) {
        label_i <- 'LFCs contrast 1'
    }
    if (is.null(label_j)) {
        label_j <- 'LFCs contrast 2'
    }

    # join results into dataframe
    cols.sub <- c('log2FoldChange', 'padj', name.col)
    df <- merge(as.data.frame(res_i)[cols.sub],
                as.data.frame(res_j)[cols.sub],
                by= name.col)
    # add significance column
    df <- df %>%
        mutate('Significance' = case_when(
                        (padj.x <= padj.thr) & (padj.y <= padj.thr) & (log2FoldChange.x * log2FoldChange.y >= 0) ~ 'Both - same LFC sign',
                        (padj.x <= padj.thr) & (padj.y <= padj.thr) & (log2FoldChange.x * log2FoldChange.y < 0) ~ 'Both - opposite LFC sign',
                        (padj.x <= padj.thr) ~ label_i,
                        (padj.y <= padj.thr) ~ label_j,
                        TRUE ~ 'None'))

    # if return.values, return the dataframe now, no need to generate the plot
    if (return.values == TRUE) {
        return(df)
    }

    # Significance as factor, to reorder in the graph
    df[['Significance']] <- factor(df[['Significance']], levels=c('None', label_j, label_i, 'Both - opposite LFC sign', 'Both - same LFC sign'))

    names(color.palette) <- c('Both - same LFC sign', 'Both - opposite LFC sign', 'None', label_i, label_j)

    p <- ggplot(df %>% arrange(Significance), aes_string(x='log2FoldChange.x', y='log2FoldChange.y',
                               color='Significance', label=name.col)) +
            geom_point(size=1) +
            theme_bw() +
            scale_color_manual(values=color.palette) +
            geom_abline(color="#333333", linetype="dashed", size=0.5, alpha=0.7) +
            geom_hline(yintercept=0, color="#333333", linetype="dashed", size=0.5, alpha=0.7) +
            geom_vline(xintercept=0, color="#333333", linetype="dashed", size=0.5, alpha=0.7) +
            xlab(label_i) +
            ylab(label_j)

    return(p)
}
