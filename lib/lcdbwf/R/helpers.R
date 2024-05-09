#' Load config file
#'
#' @param filename YAML filename to load
#'
#' @return Nested list object containing configuration
#' 
#' @details
#' For testing purposes, if the environment variable LCDB_WF_TEST exists then
#' assume we are running a test and set config$toggle$test to TRUE. This will
#' override any settings in the config YAML. Otherwise, use what is in the
#' config, and if not specified in the config then default to FALSE.
load_config <- function(filename){
  config <- yaml::yaml.load_file(filename)

  cores <- config$parallel$cores
  if (class(cores) == "character"){

    # It's possible that parallel$cores was configured with a string integer
    # rather than a bare integer. In that case we should consider it an
    # integer.
    int_cores <- as.integer(cores)
    if (!is.na(int_cores)){
      cores <- int_cores
    } else {
      cores <- Sys.getenv(cores, NA)
      if (is.na(cores)){
        warning(paste("Could not find env var", cores, "so defaulting to 1 core"))
        cores <- 1
      }
    }
  }
  # Reset value of cores based on what we've figured out above.
  config$parallel$cores <- cores

  # If test toggle is missing, default to FALSE...
  if (is.null(config$toggle$test)){ config$toggle$test <- FALSE }

  # ...but always override with env var
  is_test <- Sys.getenv('LCDB_WF_TEST', NA)
  if (!is.na(is_test)){
    config$toggle$test <- TRUE
  }
  return(config)
}


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
  obj_list <- purrr::map(var_names, function(x) eval(parse(text=x)))

  # replace names with versions without the pattern
  modified_names <- purrr::map(var_names, function(x) stringr::str_replace(x, pattern, ""))

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
lcdbwf_samplename <- function(x) {
    x <- x %>%
        stringr::str_remove_all('data/rnaseq_samples/') %>%
        stringr::str_remove_all('.cutadapt.bam') %>%
        stringr::str_split(stringr::fixed('/'), simplify=TRUE)
    x[,1]
}



#' Identify genes to label in MA plots and volcano plots
#'
#' Default is to rank by log2FoldChange and return the top 5 and bottom 5.
genes_to_label <- function(res, n=5, config){
  filtered <- res %>%
    as.data.frame() %>%
    dplyr::arrange(log2FoldChange) %>%
    dplyr::filter(!is.na(log2FoldChange))

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
    description <- res_list[[name]]$label
    tbl[[name]] <- data.frame(fn=fn, link=link, description=description)
  }
  df <- do.call(rbind, tbl) %>% dplyr::select(link, description)
  return(df)
}

#' Export Excel file of results
#'
#' One contrast per worksheet, with normalized counts and filters enabled
#'
#' @param res_list Results list object
#' @param dds_list List of dds objects
#' @param file Output .xlxs file
#'
exported_excel <- function(res_list, dds_list, file='results.xlsx'){
  wb <- openxlsx::createWorkbook()
  key <- lapply(names(res_list), function (x) res_list[[x]]$label)
  names(key) <- names(res_list)
  for (name in names(res_list)){
    res <- res_list[[name]]$res %>% as.data.frame()
    label <- res_list[[name]]$label

    # Add derived columns based on NA
    res$low_count_filtered <- (!is.na(res$padj) && is.na(res$padj))
    res$outlier <- (res$baseMean  >0 && is.na(res$padj))

    # add the normalized counts as additional columns
    dds <- dds_list[[res_list[[name]]$dds]]
    cnts <- DESeq2::counts(dds, normalized=TRUE) %>% as.data.frame()
    colnames(cnts) <- paste('normcounts_', colnames(cnts))
    cnts$gene <- rownames(cnts)
    data <- dplyr::full_join(res, cnts, by='gene')

    # Note that we can't use `label` here because there are restrictions on the
    # size and content of sheet names.
    openxlsx::addWorksheet(wb, name)
    openxlsx::writeData(wb, name, x=data, withFilter=TRUE)
  }
  openxlsx::saveWorkbook(wb, file=file, overwrite=TRUE)
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
    cnts <- DESeq2::counts(dds, normalized = TRUE, replaced = FALSE)

    # keep track of original samplenames (as we're about to add another column)
    samples <- colnames(cnts)

    # add 0.5 like plotCounts to plot 0 on log scale
    cnts <- cnts + pc

    cnts <- as.data.frame(cnts) %>% dplyr::mutate(gene=rownames(.))

    # merge with res.i and colData
    df <- dplyr::inner_join(as_tibble(cnts), tibble::as_tibble(res))

    # add label for plotting
    df <- df %>% dplyr::mutate(label = paste(gene, !!!rlang::syms(label), sep=' | '))
    # subset to sel.genes if not NULL (then keep all)
    if (!is.null(sel.genes)) {
        df <- df %>%
            dplyr::filter(gene %in% sel.genes)
    }
    # add rank
    df <- df %>%
        dplyr::mutate(rank=rank(!!sym(rank.col), ties.method='first', na.last='keep')) %>%
        tidyr::pivot_longer(all_of(samples), names_to='samplename', values_to='normalized_counts')

    # add colData, but without requiring the first column of colData to be
    # called exactly "samplename" -- instead we make a new column called
    # join.samplename that will be used just for joining (it becomes
    # "samplename" in the final joined df)

    dat <- colData(dds) %>% as.data.frame %>% dplyr::mutate(join.samplename=rownames(.))
    df <- left_join(df, dat, by=c('samplename'='join.samplename')) %>%
        tibble::as_tibble
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
my_summary <- function(res, dds, name, ...){
  dds_label <- dds
  if (class(dds) != 'character'){
    stop("expecting dds to be a string")
  }
  dds <- lcdbwf:::get_dds(dds_label)
  alpha <- metadata(res)$alpha
  lfc.thresh <- metadata(res)$lfcThreshold
  lfc.thresh <- ifelse(is.null(lfc.thresh), 0, lfc.thresh)
  notallzero <- sum(res$baseMean > 0)
  up <- sum(res$padj < alpha & res$log2FoldChange > lfc.thresh, na.rm=TRUE)
  down <- sum(res$padj < alpha & res$log2FoldChange < -lfc.thresh, na.rm=TRUE)
  filt <- sum(!is.na(res$pvalue) & is.na(res$padj))
  outlier <- sum(res$baseMean > 0 & is.na(res$pvalue))

  test <- mcols(res)['log2FoldChange', 'description']
  match <- stringr::str_match(test, ": (.*)")
  if (length(match) != 2){
    stop(paste("Expected a 1x2 matrix for matching pattern, got", match))
  }
  test <- match[1, 2]
  df <- data.frame(
    name=name,
    up=up,
    down=down,
    nonzero.vs.total=paste0(notallzero, '/', nrow(res)), 
    alpha=alpha,
    lfcThreshold=lfc.thresh,
    outliers=outlier,
    low.counts=filt,
    # adjust width.cutoff here because newline insertion causes this to return
    # a df with multiple duplicate rows
    dds=dds_label,
    design=deparse(design(dds), width.cutoff=500L),
    test=test
  )
  return(df)
}


#' Combine everything in the results list into a single table
#'
#' @param res.list Named list of lists, where each sublist contains the
#'    following names: c('res', 'dds', 'label'). "res" is a DESeqResults
#'    object, "dds" is the indexing label into the dds_list object, and "label"
#'    is a nicer-looking label to use. NOTE: backwards compatibility with older
#'    versions of lcdb-wf depends on no dds.list object being passed.
#'
#' @return data.frame
summarize_res_list <- function(res_list){
  slist <- list()
  for (name in names(res_list)){
    x <- my_summary(res_list[[name]][['res']], res_list[[name]][['dds']], name=name)
    rownames(x) <- res_list[[name]][['label']]
    slist[[name]] <- x
  }
  slist <- do.call(rbind, slist)
  rownames(slist) <- as.character(lapply(res_list, function (x) x[['label']]))
  return(slist)
}


#' Return up/down/changed genes
#'
#' @param x DESeq2 results object, or data.frame created from one
#' @param direction Direction in 'up', 'dn', 'down', 'ch', 'changed'
#' @param alpha FDR lower than this will be considered significant
#' @param thresh Log2 fold change threshold. If e.g. 2, will return < -2
#'   and/or > 2, depending on the value of "direction"
#' @param return_type If `return_type` is "rownames", return the rownames of
#'   the dataframe of the significant genes. If "results", return a subset of
#'   the original object. If "bool", then return an index of TRUE/FALSE.
#'
#' @return Character vector of rownames, bool vector, or subset of original
#'   results object depending on the value of `return_type`.
get_sig <- function(x, direction='up', alpha=0.1, lfc_thresh=0, return_type="rownames"){
  if (direction == 'up'){
    idx <- (x$padj < alpha) & (x$log2FoldChange > lfc_thresh) & (!is.na(x$padj))
  } else if (direction %in% c('down', 'dn')){
    idx <- (x$padj < alpha) & (x$log2FoldChange < -lfc_thresh) & (!is.na(x$padj))
  } else if (direction %in% c('changed', 'ch')){
    idx <- (x$padj < alpha) & (abs(x$log2FoldChange) > lfc_thresh) & (!is.na(x$padj))
  }
  if (return_type == "rownames"){
    return(rownames(x)[idx])
  } else if (return_type == "bool"){
    return(idx)
  } else if (return_type == "results"){
    return(x[idx,])
  } else {
    stop(print0("Don't know how to handle return_type='", return_type, "'."))
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


#' Compose an object to be used in downstream tools.
#'
#' @param res_list List of results objects and associated metadata. See details
#'   for format.
#' @param dds_list List of dds objects used in results
#' @param rld_list List of normalized dds objects
#' @param enrich_list List of enrichment results objects. See details for format.
#' @param degpatterns_list List of degpatterns objects
#' @param all_dds Single dds object containing all samples
#' @param all_rld Single normalized dds object containing all samples
#' @param rds_file RDS file containing lcdb-wf object. Can be used to incrementally
#'        add elements to a pre-existing run or 'sanitize' an object from a previous run.
#'        - Ignored if res_list & dds_list are specified.
#'        - If rld_list, enrich_list or degpatterns_list are also provided, these will
#'          be used to replace corresponding elements in the RDS file.
#' @param workers Number of cores to run GeneTonic conversion on
#'
#' @details
#'
#' res_list and dds_list *or* rds_file are required. `res_list` has the following format.
#'
#'    list(
#'      ko.vs.wt=list(             # names of the list are short keys
#'        res=DESeqResults object, # standard results object, may have columns for SYMBOL, etc
#'        dds="main",              # string key into dds_list
#'        label="KO vs WT"),       # Description of contrast
#'      het.vs.wt=list(
#'         ....
#'      ),
#'      ...
#'    )
#'
#' Note that values of "dds" in the res_list's items must be names in
#' `dds_list`. In this example, there is only one, "main". `dds_list` has the
#' following format:
#'
#'   list(
#'     main=DESeqDataSet,
#'     ...
#'   )
#'
#' `enrich_list` is optional. If provided, has names corresponding to
#' results names available in `res_list`. Alternatively, can have a 'res' key at the
#' second-level containing a result name available in `res_list`.
#'
#' In this example, the names are "ko.vs.wt", "het.vs.wt" & "het.vs.wt_v2". The latter two
#' both map to the `res_list` element "het.vs.wt".
#'
#' It has the following format:
#'
#'   list(
#'     ko.vs.wt=list(
#'       up=list(
#'         BP=enrichResults object,
#'         CC=enrichResults object,
#'         MF=enrichResults object,
#'         ...
#'       ),
#'       down=list(
#'         BP=enrichResults object,
#'         ...
#'       ),
#'     het.vs.wt=list(
#'       up=list(...),
#'       down=list(...)
#'     ),
#'     het.vs.wt_v2=list(
#'       res='het.vs.wt',
#'       up=list(...),
#'       down=list(...)
#'     ),
#'     ...
#'  )
#'
#'
compose_results <- function(res_list=NULL,
                            dds_list=NULL,
                            rld_list=NULL,
                            enrich_list=NULL,
                            degpatterns_list=NULL,
                            all_dds=NULL,
                            all_rld=NULL,
                            rds_file=NULL,
                            workers=1){

  if(is.null(res_list) & is.null(dds_list) & is.null(rds_file)){
    stop('Either "res_list" & "dds_list" or "rds_file" must be specified')
  } else if(is.null(res_list) | is.null(dds_list)){
    message(paste('Loading objects from RDS file:', rds_file))

    if(!file.exists(rds_file)){
      stop(paste('RDS file does not exist:', rds_file))
    }

    # get res_list & dds_list (and any others) from RDS file
    tmp <- readRDS(rds_file)

    if(!any(c('res_list', 'dds_list') %in% names(tmp))){
      stop('Object must contain "res_list" & "dds_list" elements!')
    }

    res_list <- tmp$res_list
    dds_list <- tmp$dds_list

    # plug in optional slots unless specified already
    if('rld_list' %in% names(tmp) & !is.null(rld_list)) rld_list <- tmp$rld_list
    if('enrich_list' %in% names(tmp) & !is.null(enrich_list)) enrich_list <- tmp$enrich_list
    if('degpatterns_list' %in% names(tmp) & !is.null(degpatterns_list)) degpatterns_list <- tmp$degpatterns_list
    if('all_dds' %in% names(tmp) & !is.null(all_dds)) all_dds <- tmp$all_dds
    if('all_rld' %in% names(tmp) & !is.null(all_rld)) all_rld <- tmp$all_rld
  }

  message('\n1. Processing res_list & dds_list')
  # Much of this function is just checking that the names all line up.
  res_dds_names <- unlist(lapply(res_list, function (x) x$dds))
  names(res_dds_names) <- NULL
  dds_dds_names <- names(dds_list)
  res_not_dds <- setdiff(res_dds_names, dds_dds_names)
  dds_not_res <- setdiff(dds_dds_names, res_dds_names)
  if (length(res_not_dds) > 0){
    stop(paste("\t- The following dds names are in res_list but are not found in dds_list:",
               paste(res_not_dds, collapse=', '), '\n'))
  }

  # NOTE: drop unused dds_list & rld_list objects
  if (length(dds_not_res) > 0){
    message("\t- The following dds names are in dds_list but not in res_list. These will be skipped:")
    message(paste0('\t\t', paste(dds_not_res, collapse='\n\t\t')))
    dds_list <- dds_list[ setdiff(dds_dds_names, dds_not_res) ]
    if(!is.null(rld_list)){
      rld_list <- rld_list[ setdiff(names(rld_list), dds_not_res) ]
    }
  }

  # check if rld_list was specified, if not make it
  if(is.null(rld_list)){
    message("\t- rld_list was not specified. Generating it")
    rld_list <- lapply(dds_list,
                  function(x) varianceStabilizingTransformation(x, blind=TRUE)
                )
  }

  # sanitize res_list, dds_list & rld_list
  obj <- sanitize_res_dds(res_list=res_list,
                          dds_list=dds_list,
                          rld_list=rld_list)

  # if all_dds not specified, but dds_list has length 1,
  # then use dds_list[[ 1 ]] as all_dds
  if(is.null(all_dds) & length(obj$dds) == 1){
    message('\t- all_dds was not specified, but dds_list has only 1 object. Using that instead')
    all_dds <- obj$dds[[ 1 ]]

    # if specifying all_dds, compute all_rld even if specified
    all_rld <- varianceStabilizingTransformation(all_dds, blind=TRUE)
  }

  # if all_dds is specified, but all_rld is not, compute it
  if(is.null(all_rld) & !is.null(all_dds)){
    message('\t- all_dds was specified, but not all_rld. Generating it')
    all_rld <- varianceStabilizingTransformation(all_dds, blind=TRUE)
  }

  message('\t- Generating symbol -> gene mapping from all res_list objects')
  # build symbol -> gene mapping from all res_list objects
  gene2symbol <- NULL
  gene2symbol_names <- NULL
  for(name in names(obj$res)){
    res <- obj$res[[name]]

    sidx <- which(tolower(colnames(res)) %in% 'symbol')
    gidx <- which(tolower(colnames(res)) %in% 'gene')

    gene2symbol <- c(gene2symbol, unname(res[, sidx]))
    gene2symbol_names <- c(gene2symbol_names, res[, gidx])
  }

  # remove duplicates
  idx <- !duplicated(gene2symbol)
  gene2symbol <- gene2symbol[!idx]
  names(gene2symbol) <- gene2symbol_names[!idx]

  # remove NAs
  gene2symbol[is.na(gene2symbol)] <- names(gene2symbol)[is.na(gene2symbol)]

  # replace rownames of dds_list, rld_list, all_dds, all_rld with symbol
  for(name in names(obj$dds)){
    rownames(obj$dds[[ name ]]) <- gene2symbol[ rownames(obj$dds[[ name ]]) ]
    rownames(obj$rld[[ name ]]) <- gene2symbol[ rownames(obj$rld[[ name ]]) ]
  }
  if(!is.null(all_dds)) rownames(all_dds) <- gene2symbol[ rownames(all_dds) ]
  if(!is.null(all_rld)) rownames(all_rld) <- gene2symbol[ rownames(all_rld) ]

  # plug into object
  obj[[ 'all_dds' ]] <- all_dds
  obj[[ 'all_rld' ]] <- all_rld

  if (!is.null(enrich_list)){
    message('\n2. Processing enrich_list')

    message('\t- Checking enrich_list names against res_list names')
    res_names <- names(obj$res)
    enrich_names <- names(enrich_list)
    enrich_not_res <- setdiff(enrich_names, res_names)
    if (length(enrich_not_res) > 0){
      # - if FE object name is missing in res_list, check for 'res' keys
      #   and make sure all 'res' keys are there in res_list
      # - if no 'res' keys, give error and stop
      no_res_key <- NULL
      no_res_list <- NULL
      for(name in enrich_not_res){
        if(!'res' %in% names(enrich_list[[ name ]])){
          no_res_key <- c(no_res_key, name)
        } else if(!enrich_list[['res']] %in% res_names){
          no_res_key <- c(no_res_key, name)
        }
      }

      if(length(no_res_key) > 0){
        stop(paste0("The following names are in enrich_list but do not map to res_list:\n\t", paste(no_res_key, collapse='\n\t')))
      }
    }

    # - save enrichResult objects as data.frame
    # - generate GeneTonic object
    message('\t- Converting enrich_list to genetonic objects')

    message('\t\t- Flattening enrich_list')
    # flatten enrich_list
    elem_names <- NULL
    sep <- '*'
    res_keys <- list()
    for(x in names(enrich_list)){
      for(y in names(enrich_list[[x]])){
        # NOTE: if key is 'res' save & skip
        if(y == 'res'){
          res_keys[[ x ]] <- enrich_list[[ x ]][[ 'res' ]]
          next
        }
        for(z in names(enrich_list[[x]][[y]])){
          elem_names <- c(elem_names, paste(x, y, z, sep=sep))
        }
      }
    }
    names(elem_names) <- elem_names

    message(paste('\t\t- Running conversion using', workers, 'worker(s)'))
    # run conversion
    # TODO: add check for cores if on biowulf
    flat_obj <- BiocParallel::bplapply(elem_names, function(x){
                  toks <- strsplit(x, split=sep, fixed=TRUE)[[1]]
                  if(toks[1] %in% res_names){
                    res <- obj$res[[ toks[1] ]]
                  } else if(toks[1] %in% names(res_keys)){
                    res <- obj$res[[ res_keys[[ toks[1] ]] ]]
                  } else {
                    return(NULL)
                  }

                  eres <- enrich_list[[ toks[1] ]][[ toks[2] ]][[ toks[3] ]]

                  df <- enrich_to_genetonic(eres, res)

                  df
                }, BPPARAM=BiocParallel::MulticoreParam(workers))

    message('\t- Reconstituting nested list & saving enrich_list as data frames')
    # reconstitute & clean up
    enrich_list_slim <- list()
    genetonic <- list()
    for(x in names(enrich_list)){
      enrich_list_slim[[ x ]] <- list()
      genetonic[[ x ]] <- list()

      for(y in names(enrich_list[[x]])){
        # NOTE: if key is 'res' plug in res_key & skip
        if(y == 'res'){
          enrich_list_slim[[ x ]][[ 'res' ]] <- res_keys[[ x ]]
          next
        } else {
          enrich_list_slim[[ x ]][[ y ]] <- list()
          genetonic[[ x ]][[ y ]] <- list()
        }

        for(z in names(enrich_list[[ x ]][[ y ]])){
          key <- paste(x, y, z, sep=sep)
          enrich_list_slim[[ x ]][[ y ]][[ z ]] <- enrich_list[[ x ]][[ y ]][[ z ]]@result
          genetonic[[ x ]][[ y ]][[ z ]] <- flat_obj[[key]]
        }
      }
    }

    obj[['enrich']] <- enrich_list_slim
    obj[['genetonic']] <- genetonic
  }

  if(!is.null(degpatterns_list)){
    message('\n3. Processing degpatterns_list')

    message('\t- Only keeping "normalized" slot & adding "symbol" column')
    # only keep 'normalized' slot from degpatterns object
    # & add 'symbol' column
    obj[['degpatterns']] <- lapply(degpatterns_list, function(x){
                              df <- x$normalized
                              if(!'symbol' %in% colnames(df)){
                                df$symbol <- gene2symbol[ df$genes ]
                              }
                              df
                            })
  }

  message('\nDone!')

  return(obj)
}

#' Convert enrichResult/gseaResult to GeneTonic object
#'
#' This function takes an enrichResult object and
#' DE analysis results and creates a GeneTonic object.
#'
#' @param enrich enrichResult object
#' @param res data frame with DE analysis results
#'
#' @export
enrich_to_genetonic <- function(enrich, res){
    suppressMessages({
      if(class(enrich) == 'enrichResult')
        l_gs <- shake_enrichResult(enrich)
      else if(class(enrich) == 'gseaResult')
        l_gs <- shake_gsenrichResult(enrich)
    })

    if(!'gene' %in% colnames(res)){
      if(!is.null(rownames(res))){
        res$gene <- rownames(res)
        res <- as.data.frame(res) %>% relocate(gene)
      } else {
        stop('Cannot find gene column in result data frame!')
      }
    }
    idx <- match(c('gene','symbol'), tolower(colnames(res)))
    if(length(idx) != 2){
      stop('Columns of DE results must contain "gene" & "symbol"')
    }
    anno_df <- res[,idx]
    colnames(anno_df) <- c('gene_id', 'gene_name')

    l_gs <- get_aggrscores(l_gs, res, anno_df)
    return(list(l_gs=l_gs, anno_df=anno_df))
}

#' Sanitize res_list, dds_list & rld_list for use with downstream tools
#'
#' This function makes various validation checks and sanitizes the object:
#'
#' - res_list must be a named list
#' - rownames of res_list objects cannot be NULL
#' - colData of objects cannot contain reserved column names
#' - res_list objects should have exactly one 'gene' & 'symbol' (case-insensitive).
#'   If missing, these are replaced by rownames. NA 'symbol' values
#'   are replaced by corresponding values from 'gene' column or rownames.
#' - dds_list and rld_list names must match exactly
#' - colData of dds_list & rld_list objects cannot contain reserved columns
#' - Adds a 'sample' column to colData of dds_list & rld_list objects
#' - Builds a dds_mapping object that maps res_list objects to dds_list objects
#' - Flattens res_list object which is replaced by two slots corresponding to
#'   'res' & 'label' elements.
#'
#' @param res_list List of DESeqResults objects
#' @param dds_list List of dds objects
#' @param rld_list List of normalized dds objects
#' @param reserved_cols Column names reserved for internal use. colData
#'        of dds_list or rld_list objects cannot contain these columns
#'
sanitize_res_dds <- function(res_list, dds_list, rld_list,
                             reserved_cols=c('gene', 'symbol')){
  if(is.null(names(res_list))){
    stop('"res_list" must be a named list')
  }

  if(is.null(names(dds_list))){
    stop('"dds_list" must be a named list')
  }

  for(name in names(res_list)){
    res <- res_list[[ name ]]$res

    # NOTE: rownames of res_list object cannot be NULL
    if(is.null(rownames(res))){
      stop(paste('Rownames of res_list elements cannot be NULL:', name))
    }

    # NOTE: check that a single 'gene' column exists.
    gene_idx <- grep('gene', tolower(colnames(res)))
    if(length(gene_idx) > 1){
      stop(paste('res_list elements can only have 1 "gene" column:', name))
    } else if(length(gene_idx) == 0){
      # If 'gene' column not present, replace with rownames
      message(paste('res_list element is missing a "gene" column. "rownames" will be used instead:', name))
      res$gene <- rownames(res)
      gene_idx <- grep('gene', tolower(colnames(res)))
    }

    # NOTE: check that a single 'symbol' column exists.
    symbol_idx <- grep('symbol', tolower(colnames(res)))
    if(length(symbol_idx) > 1){
      stop(paste('res_list elements can only have 1 "symbol" column:', name))
    } else if(length(symbol_idx) == 0){
      # If 'symbol' column not present, replace with rownames
      message(paste('res_list element is missing a "symbol" column. "rownames" will be used instead:', name))
      res$symbol <- rownames(res)
    } else {
      # if present, check for NA's & replace with values from 'gene' column
      na_idx <- is.na(res[, symbol_idx])
      if(sum(na_idx) > 0){
        res[na_idx, symbol_idx] <- res[na_idx, gene_idx]
      }
    }

    # plug back in to res_list
    res_list[[ name ]]$res <- res
  }

  dds_names <- names(dds_list)
  rld_names <- names(rld_list)
  if(!all(dds_names %in% rld_names)){
    stop(paste('Not all dds_list elements have matching rld_list objects:',
               setdiff(dds_names, rld_names)))
  } else if(!all(rld_names %in% dds_names)){
    stop(paste('Not all rld_list elements have matching dds_list objects:',
               setdiff(rld_names, dds_names)))
  }

  for(name in dds_names){
    dds <- dds_list[[ name ]]
    rld <- rld_list[[ name ]]

    # NOTE: colData cannot contain reserved column names
    if(any(reserved_cols %in% names(colData(dds)))){
      dds_reserved <- intersect(reserved_cols, names(colData(dds)))
      stop(paste('colData of dds_list object contains reserved column names -',
                 paste0(dds_reserved, collapse=', '), ':', name))
    }

    if(any(reserved_cols %in% names(colData(rld)))){
      rld_reserved <- intersect(reserved_cols, names(colData(rld)))
      stop(paste('colData of res_list element contains reserved column names -',
                 paste0(rld_reserved, collapse=', '), ':', name))
    }

    colData(dds)$sample <- rownames(colData(dds))
    colData(rld)$sample <- rownames(colData(rld))

    dds_list[[ name ]] <- dds
    rld_list[[ name ]] <- rld
  }

  # build res_list -> dds_list mapping to plug into degpatterns
  dds_mapping <- lapply(res_list, function(x) x$dds)
  names(dds_mapping) <- names(res_list)

  # build final object
  obj <- list(
           res=lapply(res_list, function(x) x$res),
           dds=dds_list,
           rld=rld_list,
           labels=lapply(res_list, function(x) x$label),
           dds_mapping=dds_mapping)

  return(obj)
}


#' Add cluster ID columns to res_list objects
#'
#' @param clusters DegPatterns data frame with gene -> cluster mapping
#' @param res DESeq2 results object
#' @param label Cluster column prefix
#'
#' @return DESeq2 results object with added cluster column
add.cluster.id <- function(clusters, res, label){
    if(is.null(res)) return(NULL)
    else if(is.null(clusters)){
      # add NA cluster column
      res[[paste0(label, '_cluster')]] <- NA
      return(res)
    }
    # Merges the degPattern cluster IDs `cluster` with DESeqresults `res`
    # `label` will be used to create a cluster column with a unique column name 
    # returns a DESeqresults with cluster IDs
    unq <- unique(clusters[, c('genes', 'cluster')])
    names(unq)[names(unq) == 'genes'] <- 'gene'
    merged <- merge(as.data.frame(res) %>% tibble::rowid_to_column("ID"),
                    unq, by='gene', all.x=TRUE) %>%
        arrange(ID)
    res[[paste0(label, '_cluster')]] <- merged[['cluster']]
    return(res)
}


#' Get config object from global env
#'
#' @return Config object (nested list) if in global env, otherwise NULL.
get_config <- function(){
  if ("config" %in% ls(.GlobalEnv)){
    config <- get("config", .GlobalEnv)
  } else {
    config <- NULL
  }
  return(config)
}


#' Return only the list items in `dots` that are arguments for `func`
#' 
#' Thanks to
#' https://community.rstudio.com/t/dots-vs-arg-lists-for-function-forwarding/4995/2
#' for the idea.
#'
#' @param dots List of arguments, most often created via `list(...)`
#' @param func Function that will be inspected for valid arguments
#'
#' @return List of arguments in `dots` that are valid arguments for `func`.
#'
#' @details
#'
#' Example usage is:
#'
#'  do.call("functioname", lcdbwf:::match_from_dots(list(...), functionname)))
#'
match_from_dots <- function(dots, func){
  arg <- match(names(formals(func)), names(dots))
  dots[arg[!is.na(arg)]]
}

#' Print markdown enclosed by <details> HTML tags.
#'
#' @param text Markdown text to print
#' @param summary Text to use next to dropdown arrow
folded_markdown <- function(text, summary){
  mdcat("<details>")
  mdcat("<summary>", summary, "</summary>")
  mdcat(text)
  mdcat("</details>")
}
