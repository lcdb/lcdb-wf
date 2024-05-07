#' Function to run enrichment on a res_list
#'
#' @param res_list list of DESeq2 results objects
#' @param config Config object
#' @param cores Number of cores to run it on
#' @param sep Character to separate res_list names
#'
#' @return nested list of enrichResult objects
run_enricher <- function(res_list, ontology_list, config,
                         cores=1, sep='*'){
    # This function supports running in parallel which works best with a flat
    # list; however for organizational purposese we want a nested structure. So
    # we convert between the two by collapsing nested keys for flat list, and
    # splitting the collapsed keys to reconstruct the nested.

    # make sure the sep character is not in res_list names
    if(sum(grepl(sep, names(res_list), fixed=TRUE)) > 0){
        stop(
             paste0("'res_list' names must not contain '",
                    sep, "'. Try running lcdbwf::run_enricher ",
                    "with a different 'sep' parameter."))
    }

    # Collapse names of nested res_list & ontologies
    # These will be used as keys to create a flat list structure
    n <- collapse_names(res_list=res_list,
                        config=config,
                        sep=sep)

    # run enrichment on flattened res_list
    enrich_list_flat <- BiocParallel::bplapply(n,
        function(x){
            # split name into 3 fields:
            #   comparison, direction, ontology
            toks <- unlist(strsplit(x, split=sep, fixed=TRUE))
            name <- toks[1]
            direction <- toks[2]
            ont <- toks[3]

            enrich_res <- enrich_test(
              res_list[[name]],
              direction=direction,
              TERM2GENE=ontology_list[['term2gene']][[ont]],
              TERM2NAME=ontology_list[['term2name']][[ont]],
              config=config
            )
            enrich_res
        }, BPPARAM=BiocParallel::MulticoreParam(cores))

    # create nested list structure keyed by
    # comparison, direction, ontology
    enrich_list <- list()
    for(name in names(res_list)){
        enrich_list[[name]] <- list()
        for(direction in config$functional_enrichment$directions){
            enrich_list[[name]][[direction]] <- list()
            for(ont in names(config$functional_enrichment$ontologies)){
                key <- paste(name, direction, ont, sep=sep)
                enrich_list[[name]][[direction]][[ont]] <- enrich_list_flat[[key]]
            }
        }
    }
    return(enrich_list)
}

#' Function to collapse res_list names with ontologies
#'
#' Used for converting from nested form to flattened form.
#'
#' @param res_list list of DESeq2 results objects
#' @param config Config object
#' @param collapse character string to separate the res_list names
#'
#' @return vector of strings with collapsed res_list names
collapse_names <- function(res_list, config, sep='*'){
    names <- NULL
    for(comp in names(res_list)){
        for(ch in config$functional_enrichment$directions){
            for(ont in names(config$functional_enrichment$ontologies)){
                names <- c(names,
                           paste(comp, ch, ont, sep=sep))
            }
        }
    }
    names(names) <- names
    return(names)
}

#' All-in-one enrichment function.
#'
#' Designed to not require an orgdb, and instead requires dataframes of
#' term2gene and optionally term2name.
#'
#' @param res DESeq2 results object
#' @param TERM2GENE A data.frame, first column GO ID, second column gene name.
#'   It is assumed that the gene names are the same format as those in
#'   rownames(res).
#' @param TERM2NAME A data.frame, first column GO ID, second column description.
#' @param config Config object. pvalueCutoff and qvalueCutoff will be taken from here.
#' @param direction One of "up", "down", or "changed". Will use alpha and
#'   lfc_thresh from the config.
#' @param kind One of "OR" for overrepresentation or "GSEA" for gene set
#'   enrichment analysis.
#' @param ... Additional arguments are passed on to enricher() for kind="OR" or
#'   GSEA() for kind="GSEA".
#'
#' @return An enrichResults object from
enrich_test <- function(res, TERM2GENE, TERM2NAME, config, direction, kind='OR', ...){

  if (is.null(config$main$lfc_thresh)){
    lfc_thresh <- 0
  } else {
    lfc_thresh <- config$main$lfc_thresh
  }

  if (kind == "OR"){
    genes <- get_sig(
      res$res,
      alpha=config$main$alpha,
      lfc_thresh=lfc_thresh,
      direction=direction,
      return_type="rownames"
    )

    e <- clusterProfiler::enricher(
      genes,
      TERM2GENE=TERM2GENE,
      TERM2NAME=TERM2NAME,
      pvalueCutoff=config$functional_enrichment$pvalueCutoff,
      qvalueCutoff=config$functional_enrichment$qvalueCutoff,
      ...
    )
  } else if (kind == "GSEA"){
    genes <- res$res$log2FoldChange
    names(genes) <- rownames(res$res)
    genes <- genes[!is.na(genes)]
    genes <- genes[order(genes, decreasing=TRUE)]

    e <- clusterProfiler::GSEA(
      genes,
      TERM2GENE=TERM2GENE,
      TERM2NAME=TERM2NAME,
      pvalueCutoff=config$functional_enrichment$pvalueCutoff,
      ...
    )

  } else {
    stop(paste0("Don't know how to handle enrichment type '", kind, "'."))
  }
  if (is.null(e)){
    return(e)
  }

  # portions of this copied and simplified from clusterProfiler::setReadable(),
  # but written here to avoid needing an orgdb. Also, comments added for
  # clarity.

  # keys are term IDs, values are lists of gene IDs.
  gc <- clusterProfiler::geneInCategory(e)

  # create a lookup of keytype -> label_column
  gn <- res$res[[config$annotation$label_column]]
  names(gn) <- rownames(res$res)

  # rebuild geneID column after ID conversion
  geneID <- lapply(gc, function(x){
                       paste0(gn[x], collapse='/')
            })

  # plug back into enrichResult
  eres <- e@result

  if (is(e, "gseaResult")){
    eres$core_enrichment <- unlist(geneID)
  } else {
    eres$geneID <- unlist(geneID)
  }
  e@gene2Symbol <- gn
  e@keytype <- config$annotation$label_column
  e@readable <- TRUE
  e@result <- eres
  return(e)
}


#' Get the MSigDB data for the organism provided in the config.
get_msigdb_df <- function(config){
  x <- msigdbr::msigdbr(config$annotation$genus_species)
  return(x)
}


#' Return a list of MSigDB gene sets, one per subcategory
#'
#' Note that the names of the list are concatenated gs_cat and gs_subcat, so
#' that there is a single unique key for each subcategory. We can't use just
#' subcategory because some categories (like "H") don't have subcategories.
#'
#' @param msigdb_df data.frame of MSigDB for the species, likely from
#'   get_msigdb_df() or msigdbr::msigdbr().
#'
#' @return List of dataframes, item per unique category+subcategory
#'   combination.
get_msigdb_term2gene_list <- function(msigdb_df){
  x <- msigdb_df %>%
    dplyr::mutate(geneset_key=paste(gs_cat, gs_subcat, sep="_") %>% stringr::str_replace("_$", "")) %>%
    dplyr::group_by(geneset_key)
  keys <- x %>% dplyr::group_keys()
  pieces <- x %>% dplyr::group_split()
  names(pieces) <- keys$geneset_key
  pieces <- lapply(
    pieces, function(x) {
      x %>%
        dplyr::distinct(gs_name, ensembl_gene) %>%
        as.data.frame
    }
  )
  return(pieces)
}

#' TERM2NAME for MSigDB
#'
#' @param msigdb_df data.frame of MSigDB for the species, likely from
#'   get_msigdb_df() or msigdbr::msigdbr().
get_msigdb_term2name <- function(msigdb_df){
  x <- msigdb_df %>% dplyr::distinct(gs_name, gs_description) %>% as.data.frame
  return(x)
}



#' Alternative implementation of get_go_term2gene, stored here for historical
#' comparison purposes.
get_go_term2gene_alt <- function(orgdb, keytype){
  kk <- keys(orgdb, keytype=keytype)
  goAnno <- AnnotationDbi::select(orgdb, keys=kk, keytype=keytype, columns=c("GOALL", "ONTOLOGYALL"))
  goAnno <- unique(goAnno[!is.na(goAnno[,1]), ])
  lst <- list(
      MF=goAnno %>% dplyr::filter(ONTOLOGYALL=="MF") %>% dplyr::select(GOALL, !!keytype),
      CC=goAnno %>% dplyr::filter(ONTOLOGYALL=="CC") %>% dplyr::select(GOALL, !!keytype),
      BP=goAnno %>% dplyr::filter(ONTOLOGYALL=="BP") %>% dplyr::select(GOALL, !!keytype)
  )
  return(lst)

}

#' Get the species name that works with the KEGG database
#'
#' The KEGG pathway needs the species name in a
#' specific format, e.g. for 'Homo sapiens' the
#' KEGG version would be 'hsa'.
#'
#' @param config Config object
#'
#' @return KEGG-compatible species name
get_kegg_species <- function(config){
    species <- config$annotation$genus_species
    species_split <- unlist(strsplit(species, "\\s+"))
    kegg_species <- paste0(tolower(substr(species_split[1],1,1)),
                           substr(species_split[2],1,2))
    return(kegg_species)
}

#' All-in-one function to get ontology information
#'
#' @param config Config object.
#'
#' @return List of lists where top-level elements are ontology
#'      names. Each element is a list with TERM2GENE & TERM2NAME
#'      data frames
get_ontology_list <- function(config){

    orgdb <- lcdbwf:::get_orgdb(config)
    keytype <- config$annotation$keytype

    # This is the method used by clusterProfiler internally. See
    # get_go_term2gene_alt for a different implementation.
    goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
    go2gene <- suppressMessages(
      AnnotationDbi::mapIds(
        orgdb, keys=names(goterms), column=keytype, keytype="GOALL",
        multiVals='list')
    )
    goAnno <- stack(go2gene)
    colnames(goAnno) <- c(keytype, "GOALL")
    goAnno <- unique(goAnno[!is.na(goAnno[,1]), ])
    goAnno$ONTOLOGYALL <- goterms[goAnno$GOALL]

    # get GO descriptions
    go2name <- AnnotationDbi::select(GO.db::GO.db,
                    keys=keys(GO.db::GO.db, "GOID"),
                    c("GOID", "TERM"))

    # Split up the dataframe into a list, one per annotation.
    go2gene <- list(
        MF=goAnno %>% dplyr::filter(ONTOLOGYALL=="MF") %>% dplyr::select(GOALL, !!keytype),
        CC=goAnno %>% dplyr::filter(ONTOLOGYALL=="CC") %>% dplyr::select(GOALL, !!keytype),
        BP=goAnno %>% dplyr::filter(ONTOLOGYALL=="BP") %>% dplyr::select(GOALL, !!keytype)
    )

    # We need to assign each key to its respective term2name dataframe (or NULL if
    # none)
    ontology_list <- list(term2gene=go2gene,
                          term2name=lapply(go2gene,
                            function(x) go2name %>% dplyr::filter(GOID %in% x$GOALL)))


    # download KEGG information
    kegg_species <- lcdbwf:::get_kegg_species(config)

    # NOTE: KEGG API changes breaks clusterProfiler enrichKEGG
    #       and download_KEGG functions pre v4.7.2. Using
    #       internal patched versions in its place here.
    #       Can switch back to clusterProfiler version later
    #       if needed, so leaving the alternate command here.
    #
    #kegg_list <- clusterProfiler::download_KEGG(kegg_species,
    #                            keyType='kegg')
    kegg_list <- get_KEGG_info(kegg_species)

    term2gene <- kegg_list[['term2gene']]
    term2name <- kegg_list[['term2name']]

    # get pathway to gene ID mapping
    term2id <- NULL
    while(is.null(term2id)){
      term2id <- tryCatch(
        suppressMessages(
          AnnotationDbi::mapIds(
              orgdb, keys=term2gene[,2],
              column=keytype, keytype="ENTREZID",
              multiVals='first')
        ),
        error = function(e){ NULL }
      )

      # if null, convert IDs to 'ncbi-geneid' ('ENTREZID') first
      if(is.null(term2id)){
        idconv <- download_KEGG_db('ncbi-geneid', kegg_species,
                                   'conv')
        idconv[, 1] <- gsub('.+\\:', '', idconv[,1])
        idconv[, 2] <- gsub('.+\\:', '', idconv[,2])
        rownames(idconv) <- idconv[,1]

        term2gene[, 2] <- idconv[term2gene[,2], 2]
      }
    }

    term2gene[,2] <- term2id[term2gene[,2]]
    colnames(term2gene) <- c('KEGG', keytype)
    colnames(term2name) <- c('KEGG', 'Description')

    # add kegg info to ontology_list
    ontology_list$term2gene[['KEGG']] <- term2gene
    ontology_list$term2name[['KEGG']] <- term2name

    # Get all of MSigDB, although we may only use subsets of it.
    # This can take up a lot of memory on CI/CD, so we only do this if not doing
    # a test.
    if (!config$toggle$test){
      msigdb_df <- get_msigdb_df(config)
      msigdb_term2gene_list <- get_msigdb_term2gene_list(msigdb_df)
      ontology_list$term2gene <- c(ontology_list$term2gene,
                                   msigdb_term2gene_list)

      # MSigDB term names are very long and don't convey
      # that much more information than the names. So,
      # setting that to NULL here
      ontology_list$term2name <- c(ontology_list$term2name,
                                   lapply(msigdb_term2gene_list,
                                          function(x) NULL))
    }

    return(ontology_list)
}

#' Download KEGG info from API
#'
#' This function is a simplified version combining clusterProfiler
#' functions 'kegg_rest' & 'kegg_link' to reflect KEGG API
#' changes.
#'
#' @param target_db KEGG species db name, e.g. 'hsa' for human
#' @param source_db KEGG source db name, e.g. 'pathway'.
#' @param type type of information to download. Can be 'term2gene', 'term2name' or 'conv'
#'
download_KEGG_db <- function(target_db, source_db, type){
  if(type == 'term2gene') type <- 'link'
  else if(type == 'term2name') type <- 'list'
  else if(type == 'conv') type <- 'conv'
  else {
    stop('Unrecognized "type": must be "term2gene", "term2name" or "conv"')
  }

  # first get: pathway -> gene id mapping
  url <- paste("https://rest.kegg.jp", type,
                target_db, source_db, sep = "/")

  # download data to tempfile and save to data frame
  f <- tempfile()
  dl <- tryCatch(
          utils::download.file(url, quiet=TRUE,
                               method='libcurl',
                               destfile=f),
                 error = function(e) NULL
        )

  if (is.null(dl)) {
      message("Failed to download KEGG data.")
      return(NULL)
  }
  content <- readLines(f)
  content %<>% strsplit(., "\t") %>% do.call("rbind", .)
  res <- data.frame(from = content[, 1], to = content[, 2])

  # next get pathway -> pathway name mapping

  return(res)
}

#' Build annotation list from KEGG info
#'
#' This wrapper function downloads term -> gene & term -> name
#' mappings from KEGG and returns as a list of data frames.
#'
#' @param species KEGG species name, e.g. 'hsa'
#'
get_KEGG_info <- function(species){
  term2gene <- download_KEGG_db(species, 'pathway', 'term2gene')
  term2name <- download_KEGG_db('pathway', species, 'term2name')

  # clean up term2gene prefixes
  term2gene[, 'from'] <- gsub('.+\\:', '', term2gene[, 'from'])
  term2gene[, 'to'] <- gsub('.+\\:', '', term2gene[, 'to'])

  return(
    list(term2gene=term2gene,
         term2name=term2name)
  )
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
    filename.orig <- file.path(cprof.folder, paste0(label, '.tsv'))
    write.table(res, file=filename.orig, sep='\t', quote=FALSE, row.names=FALSE)
    filename.split <- file.path(cprof.folder, paste0(label, '_split.tsv'))
    res.split <- split.clusterProfiler.results(res)
    write.table(res.split, file=filename.split, sep='\t', quote=FALSE, row.names=FALSE)
    return(list(orig=filename.orig, split=filename.split))
}

#' Apply a function to a nested enrichResults list
#'
#' @param enrich_list Nested list. First level keys are results names. Second
#'   level keys are directions (up, down, changed). Third level keys are
#'   ontology labels ("ont").
#' @param func Function to apply to every identified enrichResult object
#' @param send_names If TRUE, `func` will also be provided with arguments for
#'   name, direction, and ont. This can be useful if the function needs to
#'   behave differently depending on the object.
#'
#' @return A list of the same shape (same nesting, same keys) but with the
#' "leaves" replaced with whatever `func` returns.
#' 
#' @details
#' The results from many functional enrichment runs across many contrasts can
#' be tedious to work with. This function can be used to apply a function to
#' each "leaf" enrichResults object in the list.
#'
#' The list is expected to have the structure:
#'
#'   list(
#'     contrast_name=list(
#'       direction=list(
#'         ontology_name=enrichResults)
#'       )
#'   )
#'
#' That is, `enrich_list[["ko.vs.wt"]][["up"]][["BP"]]` will get the GO
#' Biological Process enrichment results from the upregulated genes of contrast
#' "ko.vs.wt".
#'
enrich_list_lapply <- function(enrich_list, func, send_names=FALSE, ...){
  out <- list()
  for (name in names(enrich_list)){
    out[[name]] <- list()
    for (direction in names(enrich_list[[name]])){
      out[[name]][[direction]] <- list()
      for (ont in names(enrich_list[[name]][[direction]])){
        x <- enrich_list[[name]][[direction]][[ont]]
        if (send_names){
          out[[name]][[direction]][[ont]] <- func(x, name=name, direction=direction, ont=ont, ...)
        } else {
          out[[name]][[direction]][[ont]] <- func(x, ...)
        }
      }
    }
  }
  return(out)
}


#' Truncate the names of an enrichment results object
#'
#' @param obj DOSE::enrichResult object
#' @param truncate_to Max number of characters in a label
#'
#' @return enrichResult object with labels truncated
truncate <- function(obj, truncate_to){
  f <- function(x){
    if (nchar(x) > truncate_to){
      x <- paste0(substr(x, 1, truncate_to), '...')
    }
    return(x)
  }
  obj@result$Description <- sapply(obj@result$Description, f)
  obj@result$Description <- make.unique(obj@result$Description)
  return(obj)
}


#' Wrapper for enrichplot::dotplot
#'
#' @param enrich_res enrichResult object
#' @param config Config object
#' @param name, direction, ont These may be passed in when iterating over
#'   nested list of results
#' @param truncate_to Truncate labels to this long
dotplots <- function(enrich_res, config, name=NULL, direction=NULL, ont=NULL, truncate_to=50){
  if (!is.null(name)){
    title <- paste(name, direction, ont, sep=", ")
  } else {
    title <- ""
  }
  enrich_res <- truncate(enrich_res, truncate_to)
  p <- do.call(enrichplot::dotplot, c(enrich_res, config$plotting$dotplot_args)) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5, size=15, face='bold'))
  return(p)
}

#' Wrapper for enrichplot::emapplot
#'
#' @param enrich_res enrichResult object
#' @param config Config object
#' @param name, direction, ont These may be passed in when iterating over
#'   nested list of results
#' @param truncate_to Truncate labels to this long
emapplots <- function(enrich_res, config, name=NULL, direction=NULL, ont=NULL, truncate_to=50){
  if (!is.null(name)){
    title <- paste(name, direction, ont, sep=", ")
  } else {
    title <- ""
  }
  enrich_res <- truncate(enrich_res, truncate_to)
  enrich_res <- enrichplot::pairwise_termsim(enrich_res)
  p <- do.call(enrichplot::emapplot, c(enrich_res, config$plotting$emapplot_args)) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5, size=15, face='bold'))
  return(p)
}

#' Wrapper for enrichplot::cnetplot
#'
#' @param enrich_res enrichResult object
#' @param config Config object
#' @param name, direction, ont These may be passed in when iterating over
#'   nested list of results
#' @param truncate_to Truncate labels to this long
cnetplots <- function(enrich_res, config, name=NULL, direction=NULL, ont=NULL, truncate_to=50){
  if (!is.null(name)){
    title <- paste(name, direction, ont, sep=", ")
  } else {
    title <- ""
  }
  enrich_res <- truncate(enrich_res, truncate_to)
  p <- do.call(enrichplot::cnetplot, c(enrich_res, config$plotting$cnetplot_args)) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5, size=15, face='bold'))
  return(p)
}

#' List keys for MSigDB that can be added to
#' config$functional_enrichment$ontologies
#"
#' @return List of keys using the same concatenation of gs_cat and gs_subcat
#'   that is used in the get_msigdb_term2gene_list function.
available_msigdb_keys <- function(){
  df <- msigdbr::msigdbr_collections() %>%
    as.data.frame %>%
    mutate(x=paste(gs_cat, gs_subcat, sep='_') %>% str_replace('_$', ''))
  return(df$x)
}
