# ------------------------------------------------------------------------------
# Functions for working with annotations and OrgDbs
# ------------------------------------------------------------------------------

#' Get the OrgDb for the specified organism, using the cached AnnotationHub.
#'
#' @param config List containing the following named items:
#'  - annotation_genus_species: e.g., "Homo sapiens"
#'  - hub_cache: Optional path, relative to the working directory, in which
#'               to store the downloaded OrgDb database
#'  - annotation_key_override: Optional key that if specified overrides
#'    annotation_genus_species. Use this to be exact about the version of the
#'    OrgDb you want to use.
#' @return OrgDb object
get.orgdb <- function(config){

    annotation_genus_species <- config[['annotation_genus_species']]
    annotation_key_override <- config[['annotation_key_override']]
    cache <- config[['hub_cache']]

    # Workaround to allow AnnotationHub to use proxy. See
    # https://github.com/Bioconductor/AnnotationHub/issues/4, and thanks
    # Wolfgang!
    proxy <- Sys.getenv('http_proxy')
    if (proxy == ""){
        proxy <- NULL
    }

    if (!dir.exists(cache)){
        dir.create(cache, recursive=TRUE)
    }

    ah <- AnnotationHub(hub=getAnnotationHubOption('URL'),
         cache=cache,
         proxy=proxy,
         localHub=FALSE
    )

    if (is.null(annotation_key_override)) {
        ah.query <- query(ah, "OrgDb")
        ah.query.speciesmatch <- grepl(paste("^", annotation_genus_species, "$", sep=""), ah.query$species)
        ah.query.which <- which(ah.query.speciesmatch)
        stopifnot(length(ah.query.which) > 0) #require at least one match
        if (length(ah.query.which) > 1) { #warn of duplicate matches
           print("WARNING: found multiple candidate species in AnnotationHub: ");
           print(ah.query.speciesmatch)
        }
        annotation_key <- names(ah.query)[ah.query.which[1]]
    } else {
        annotation_key <- annotation_key_override
    }

    orgdb <- ah[[annotation_key]]
    return(orgdb)
}

#' Attach additional OrgDb information to results objects within a list
#'
#' This can also fill in NAs in columns specified in the config with rownames
#' (e.g., fill in missing symbols with Ensembl IDs)
#'
#' @param res_list List of DESeqResults objects
#' @param config Full config object, at least containing config$orgdb
#'
#' @return List of same results objects, but each one with additional columns
#'   attached as specified in the config
attach_extra <- function(res_list, config){
  #' Ensure everything has the same gene IDs; if so then we'll grab just the
  #' gene IDs from the first one to use as our keys.
  nms <- lapply(res_list, function(x) rownames(x$res))
  if (!all(sapply(nms, function(x) all(x==nms[[1]])))){
    stop("Gene IDs differ between results objects")
  }
  keys <- nms[[1]]

  orgdb <- lcdbwf::get.orgdb(config$orgdb)

  # Create a dataframe mapping gene IDs to the various configured columns
  lookups <- list()
  for (col in config$orgdb$orgdb_columns){
    lookups[[col]] <- mapIds(orgdb, keys=keys, column=col, keytype=config$orgdb$keytype, multiVal='first')
    if (col %in% config$orgdb$fill){
      lookups[[col]] <- ifelse(is.na(lookups[[col]]), keys, lookups[[col]])
    }
  }
  lookups <- data.frame(lookups)

  # Use that dataframe to attach additional columns to each results object
  for (name in names(res_list)){
    res <- res_list[[name]]$res
    orig_colnames <- colnames(res)
    for (col in config$orgdb$orgdb_columns){
      res[[col]] <- lookups[[col]]
    }
    res$gene <- rownames(res)
    res <- res[, c('gene', config$orgdb$orgdb_columns, orig_colnames)]
    res_list[[name]]$res <- res
  }
  return(res_list)
}
