

# All-in-one function that does functional enrichment on a single DESeqResults
# object.
functional_enrichment <- function(res, orgdb, config, ontology, kind='OR', direction, lfc_thresh=0, alpha=0.1){
  allowed_ontologies <- c('BP', 'MF', 'CC', 'KEGG', 'Reactome')
  if (!(ontology %in% allowed_ontologies)){
    stop(paste("No current support for ontology ", ontology))
  }

  # Overrepresentation tests only work on sets
  if (kind == 'OR'){
    subset_res <- get.sig(res, direction=direction, alpha=0.1, lfc.thresh=lfc_thresh, return.names=FALSE)
    genes <- rownames(subset_res)

    if (ontology %in% c('BP', 'MF', 'CC')){
      or_res <- enrichGO(
        gene=genes,
        universe=config$functional_enrichment$go_uni,
        keyType=config$functional_enrichment$go_keytype,
        OrgDb=orgdb,
        ont=ontology,
        pvalueCutoff=config$functional_enrichment$pvalueCutoff,
        qvalueCutoff=config$functional_enrichment$qvalueCutoff,
        readable=TRUE)

    } else if (ontology == 'KEGG'){

      # Note that IDs will need to be converted later when we have full access to
      # the res_list
      or_res <- enrichKEGG(
        gene=genes,
        universe=config$functional_enrichment$path_uni,
        organism=config$functional_enrichment$kegg_organism,
        keyType=config$functional_enrichment$kegg_keytype,
        pvalueCutoff=config$functional_enrichment$pvalueCutoff,
        qvalueCutoff=config$functional_enrichment$qvalueCutoff
      )

      or_res <- clusterProfiler::setReadable(or_res, OrgDb=orgdb, keyType=confi$functional_enrichment$kegg_keytype)

    } else if (ontology == 'Reactome') {
      if (!config$toggle$reactome) {
        or_res <- enrichPathway(
          gene=genes,
          universe=config$functional_enrichment$pathway_uni,
          organism=config$functional_enrichment$reactome_organism,
          pvalueCutoff=config$functional_enrichment$pvalueCutoff,
          qvalueCutoff=config$functional_enrichment$qvalueCutoff,
          readable=TRUE)
      } else {
        stop("Reactome is toggled off in the config")
      }
    }

  } else if (kind == 'GSEA') {
    genes <- rownames(res)
  }
  return(or_res)
}




#' Run overrepresentation analysis
#'
#' @param genes.df Nested list of lists, see details below.
#' @param orgdb OrgDb containing GO and KEGG terms
#' @param orgdb_config Config list (documented elsewhere)
#' @param truncate_labels_to Long labels will be truncated to this many
#'        characters to make plots more readable
#'
#' @details
#' `genes.df` is a nested list of lists. The first level represents contrasts
#' (condition1 vs condition2, etc) and within each contrast are the selections
#' (typically "up" and "dn"). The values are character vectors of gene IDs. For
#' example:
#'
#' genes.df <- list(
#'    contrast1=list(up=c('gene1', 'gene2'),
#'                   dn=c('gene8', 'gene9', 'gene10')
#'                   ),
#'    contrast2=list(...
#'    )
#' )
#'
#' `orgdb_config` is a nested list containing configuration information and is
#' documented elsewhere.
run_functional_enrichment <- function(genes.df, orgdb, orgdb_config, truncate_labels_to=50){

    # all.enrich will hold all enrichment results. It is a nested
    # list-of-lists-of-lists, with the structure:
    #
    #   all.enrich[[contrast]][[direction]][[enrichment label]]

    all.enrich <- list()
    enrich.files <- list()

    # loop over each contrast
    for(comp in names(genes.df)){

        all.enrich[[comp]] <- list()
        enrich.files[[comp]] <- list()

        # loop over selection list
        for(sel in sel.names){

            all.enrich[[comp]][[sel]] <- list()
            enrich.files[[comp]][[sel]] <- list()

            gene <- genes.df[[comp]][[sel]][[functional_enrichment_config[['go_keytype']]]]

            for (ont in c('CC', 'BP', 'MF')){
                go.res <- enrichGO(
                    gene=gene,
                    universe=functional_enrichment_config[['go_uni']],
                    keyType=functional_enrichment_config[['go_keytype']],
                    OrgDb=orgdb,
                    ont=ont,
                    pvalueCutoff=functional_enrichment_config[['pvalueCutoff']],
                    qvalueCutoff=functional_enrichment_config[['qvalueCutoff']],
                    readable=TRUE)

                res.label <- paste('GO', ont, sep='_')
                all.enrich[[comp]][[sel]][[ont]] <- go.res
                if (!is.null(go.res)){
                    enrich.files[[comp]][[sel]][[ont]] <- write.clusterprofiler.results(go.res, cprof.folder, paste0(res.label, '_', comp, '_', sel))
                }
            }

            # perform KEGG enrichment
            res.label <- 'KEGG'
            kegg.res <- enrichKEGG(
                gene=gene,
                universe=functional_enrichment_config$path_uni, organism=functional_enrichment_config$kegg_organism,
                keyType=functional_enrichment_config$kegg_keytype,
                pvalueCutoff=functional_enrichment_config$pvalueCutoff,
                qvalueCutoff=functional_enrichment_config$qvalueCutoff
            )


            if (!is.null(kegg.res)){
                # convert uniprot IDs to readable symbols. This needs to be done separately
                # in the case of KEGG
                if(functional_enrichment_config$id_convert){
                  id.vec <- genes.df[[comp]][[sel]]$SYMBOL
                  names(id.vec) <- genes.df[[comp]][[sel]]$ENTREZID
                  kegg.res.genes <- kegg.res@result$geneID
                  for(j in 1:length(kegg.res.genes)){
                    id.split <- strsplit(kegg.res.genes[j], "/")[[1]]
                    temp <- paste(id.vec[id.split], collapse="/")
                    kegg.res@result$geneID[j] <- temp
                  }
                }


                all.enrich[[comp]][[sel]][['kegg']] <- kegg.res
                enrich.files[[comp]][[sel]][['kegg']] <- write.clusterprofiler.results(kegg.res, cprof.folder, paste0(res.label, '_', comp, '_', sel))
            }

            if (toggle$run_reactome) {
                # NOTE: Install OrgDb--------------------------------------------------
                #   ReactomePA requires the relevant OrgDb package to be installed
                #   and so is disabled by default.
                res.label <- 'Reactome'
                reactome.res <-enrichPathway(
                    gene=gene,
                    universe=functional_enrichment_config[['path_uni']],
                    organism=functional_enrichment_config[['reactome_organism']],
                    pvalueCutoff=functional_enrichment_config[['pvalueCutoff']],
                    qvalueCutoff=functional_enrichment_config[['qvalueCutoff']],
                    readable=TRUE)
                all.enrich[[comp]][[sel]][['reactome']] <- reactome.res

                if (!is.null(reactome.res)){
                    enrich.files[[comp]][[sel]][['reactome']] <- write.clusterprofiler.results(reactome.res, cprof.folder, paste0(res.label, '_', comp, '_', sel))
                }
            }
        }
    }

    # Truncate names to be <50 characters. This helps make the plot labels more
    # reasonable.
    all.enrich <- lcdbwf::truncate_names(all.enrich, truncate_labels_to)

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


truncate_names <- function(all_enrich, maxcatlen=50){
    for(comp in names(all_enrich)){
        for(name in names(all_enrich[[comp]])){
            for(go in names(all_enrich[[comp]][[name]])){
                all_enrich[[comp]][[name]][[go]]@result$Description <- substr(
                    all_enrich[[comp]][[name]][[go]]@result$Description, 1, maxcatlen
                )
            }
        }
    }
    return(all_enrich)
}

