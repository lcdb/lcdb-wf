# -----------------------------------------------------------------------------
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


# RMarkdown utilities --------------------------------------------------------

#' Find named chunks that match a pattern
#'
#' Useful for depending on chunks that are named a certain way.
#'
#' @param pattern Regular expression used to find chunks.
#'
#' @return List of chunk names matching pattern
find_chunk_names <- function(pattern){
  knitr::all_labels()[grepl(pattern, knitr::all_labels())]
}

#' Find objects whose names match a pattern and store in a list
#'
#' The pattern will be used to search within the global environment. The names
#' of the list will be the variable names with the matched pattern removed.
#'
#' @param pattern Pattern to provide to grepl
#' @param fixed Passed to grepl()
#'
#' @return List of identified objects
collect_objects <- function(pattern, fixed=FALSE){

  # Without .GlobalEnv we only catch variables within this function. Can also
  # use set.frame() to specify a frame relative to this one, but .GlobalEnv is
  # probably the most generally useful
  var_names <- ls(.GlobalEnv)
  var_names <- var_names[grepl(pattern, var_names, fixed=fixed)]

  # names are variable names; values are objects
  obj_list <- map(var_names, function(x) eval(parse(text=x)))

  # replace names with versions without the pattern
  modified_names <- map(var_names, function(x) str_replace(x, pattern, ""))

  # If there was a wildcard in the pattern there is a risk that the modified
  # names are no longer unique
  if (length(unique(modified_names)) != length(var_names)){
    stop(paste("Found duplicate names after removing pattern", pattern))
  }
  names(obj_list) <- modified_names
  return(obj_list)
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



#' Identify genes to label in MA plots and volcano plots
#'
#' Default is to rank by log2FoldChange and return the top 5 and bottom 5.
genes_to_label <- function(res, n=5, config){
  filtered <- res %>%
    as.data.frame() %>%
    arrange(log2FoldChange) %>%
    filter(!is.na(log2FoldChange))

  genes <- c(
    head(filtered, n)[[config$annotation$label_column]],
    tail(filtered, n)[[config$annotation$label_column]]
  )
  return(genes)
}


#' Save TSVs to ourput directory, named after the names in the list.
#'
#' @param res_list "reslist" object
#' @param directory Directory in which to save results
exported_tsvs <- function(res_list, directory="results"){
  dir.create(directory, showWarnings=FALSE, recursive=TRUE)
  tbl <- list()
  for (name in names(res_list)){
    fn <- file.path(directory, paste0(name, '.tsv'))
    write.table(res_list[[name]]$res, file=fn, row.names=FALSE, sep='\t', quote=FALSE)
    link <- paste0('[', fn, '](', fn, ')')
    description <- paste0('results for ', res_list[[name]]$label)
    tbl[[name]] <- data.frame(fn=fn, link=link, description=description)
  }
  df <- do.call(rbind, tbl) %>% dplyr::select(link, description)
  return(df)
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

    # keep track of original samplenames (as we're about to add another column)
    samples <- colnames(cnts)

    # add 0.5 like plotCounts to plot 0 on log scale
    cnts <- cnts + pc

    cnts <- as.data.frame(cnts) %>% mutate(gene=rownames(.))

    # merge with res.i and colData
    df <- inner_join(as_tibble(cnts), as_tibble(res))

    # add label for plotting
    df <- df %>% mutate(label = paste(gene, !!!syms(label), sep=' | '))
    # subset to sel.genes if not NULL (then keep all)
    if (!is.null(sel.genes)) {
        df <- df %>%
            filter(gene %in% sel.genes)
    }
    # add rank
    df <- df %>%
        mutate(rank=rank(!!sym(rank.col), ties.method='first', na.last='keep')) %>%
        pivot_longer(all_of(samples), names_to='samplename', values_to='normalized_counts')

    # add colData, but without requiring the first column of colData to be
    # called exactly "samplename" -- instead we make a new column called
    # join.samplename that will be used just for joining (it becomes
    # "samplename" in the final joined df)

    dat <- colData(dds) %>% as.data.frame %>% mutate(join.samplename=rownames(.))
    df <- left_join(df, dat, by=c('samplename'='join.samplename')) %>%
        as_tibble
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
    rownames(slist) <- as.character(lapply(res.list, function (x) x[['label']]))
    return(slist)
}


#' Return index of up/down/changed genes
#'
#' @param x DESeq2 results object, or data.frame created from one
#' @param direction Direction in 'up', 'dn', 'down', 'ch', 'changed'
#' @param alpha FDR lower than this will be considered significant
#' @param thresh Log2 fold change threshold. If e.g. 2, will return < -2
#'   and/or > 2, depending on the value of "direction"
#' @param return.names If TRUE, returns the rownames of selected genes; if
#'   FALSE return boolean index of length(x)
#' @param column If provided and return.names is TRUE, then return this column
#' instead of the default rownames
#'
#' @return Character vector of rownames (if return.names=TRUE) or boolean vector of genes selected.
get.sig <- function(x, direction='up', alpha=0.1, lfc.thresh=0, return.names=TRUE, column=NULL){
    if (direction == 'up'){
        idx <- (x$padj < alpha) & (x$log2FoldChange > lfc.thresh) & (!is.na(x$padj))
    } else if (direction %in% c('down', 'dn')){
        idx <- (x$padj < alpha) & (x$log2FoldChange < -lfc.thresh) & (!is.na(x$padj))
    } else if (direction %in% c('changed', 'ch')){
        idx <- (x$padj < alpha) & (abs(x$log2FoldChange) > lfc.thresh) & (!is.na(x$padj))
    }
    if (return.names){
        if (!is.null(column)){
          return(x[idx, column])
        } else {
          return(rownames(x)[idx])
        }
    } else {
        return(idx)
    }
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

