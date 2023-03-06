# ------------------------------------------------------------------------------
# Functions for working with annotations and OrgDbs
# ------------------------------------------------------------------------------


#' Get the AnnotationHub object, using settings from config
#'
#' @param config Config object
#' @param localHub, force, cache Arguments passed on to AnnotationHub(); these
#'   can be used to override what's in the config.
get_annotation_hub <- function(config, localHub=NULL, force=NULL, cache=NULL){

  if (missing(localHub)) localHub <- config$annotation$localHub
  if (is.null(localHub)) localHub <- FALSE
  if (missing(force)) force <- config$annotation$force
  if (missing(cache)) cache <- config$annotation$cache
  if (is.null(cache)) cache <- AnnotationHub::getAnnotationHubOption("CACHE")

  proxy <- Sys.getenv('http_proxy')
  if (proxy == ""){
      proxy <- NULL
  }

  if (!is.null(cache)){
    if (!dir.exists(cache)) dir.create(cache, recursive=TRUE)
  }

  ah <- AnnotationHub::AnnotationHub(
    hub=getAnnotationHubOption('URL'),
    proxy=proxy,
    localHub=localHub,
    cache=cache
  )

  # AnnotationHub uses a safe permissions approach, setting the AnnotationHub
  # lock file to be only visible by the creating user and the cache database to
  # be read-only for the group. However, this can cause permission errors when
  # working in a group setting. If this setting is TRUE, then the permissions
  # will be set on BiocFileCache.sqlite and BiocFileCache.sqlite.LOCK to be
  # read/write for both user and group.
  if (config$main$group_permissions){
    files <- dir(cache, full.names=TRUE)
    files <- files[grep('BiocFileCache.sqlite', files)]
    Sys.chmod(files, mode="0660", use_umask=TRUE)
  }

  return(ah)
}


#' Convenience wrapper function for get_annotation_db()
get_orgdb <- function(...){get_annotation_db(..., dbtype="OrgDb")}


#' Convenience wrapper function for get_annotation_db()
get_txdb <- function(...){get_annotation_db(..., dbtype="TxDb")}


#' Get the OrgDb for the specified organism, using the cached AnnotationHub.
#'
#' @param config List containing the following named items:
#'  - annotation_genus_species: e.g., "Homo sapiens"
#'
#'  - hub_cache: Optional path, relative to the working directory, in which
#'    to store the downloaded OrgDb database
#'
#'  - annotation_key_override: Optional key that if specified overrides
#'    annotation_genus_species. Use this to be exact about the version of the
#'    OrgDb you want to use.
#'
#' @return OrgDb object
get_annotation_db <- function(config, dbtype, genus_species=NULL, orgdb_key_override=NULL, txdb_key_override=NULL, cache=NULL){

    if (missing(genus_species)) genus_species <- config$annotation$genus_species
    if (missing(cache)) cache <- config$annotation$hub_cache
    if (missing(orgdb_key_override)) orgdb_key_override <- config$annotation$orgdb_key_override
    if (missing(txdb_key_override)) txdb_key_override <- config$annotation$txdb_key_override

    ah <- lcdbwf:::get_annotation_hub(config, cache=cache)

    # If an override was provided, immediately return the corresponding db.
    if (!is.null(orgdb_key_override) & dbtype == "OrgDb") return(ah[[orgdb_key_override]])
    if (!is.null(txdb_key_override) & dbtype == "TxDb") return(ah[[txdb_key_override]])

    ah.query <- AnnotationHub::query(ah, c(genus_species, dbtype))


    # Ensure that the species matches exactly
    if (!all(ah.query$species == genus_species)){ 
      stop(
        "Multiple species matched '", genus_species, "' ",
        "in the AnnotationHub search. Please use a more precise name, ",
        "or consider manually searching an using the orgdb_key_override ",
        "argument."
      )
    }

    hits <- mcols(ah.query) %>%
      as.data.frame() %>%
      dplyr::arrange(desc(rdatadateadded)) %>%
      dplyr::filter(rdataclass==dbtype)

    db_to_use <- hits[1,]

    if (nrow(hits) < 1){
      warning(paste("Found multiple", dbtype, "hits for", genus_species,
                    ". Using latest from", db_to_use$rdatadateadded))
    }

    annotation_key <- rownames(db_to_use)
    return(ah[[annotation_key]])
}

#' Attach additional OrgDb information to results objects within a list
#'
#' This can also fill in NAs in columns specified in the config with rownames
#' (e.g., fill in missing symbols with Ensembl IDs). See config documentation
#' for details.
#'
#' @param res_list List of DESeqResults objects
#' @param config Full config object, at least containing config$orgdb
#' @param use_orgdb Boolean to use or bypass orgDb extra columns
#'
#' @return List of same results objects, but each one with additional columns
#'   attached as specified in the config
attach_extra <- function(res_list, config, force_intersect, use_orgdb=TRUE){

  if (missing(force_intersect)) force_intersect <- config$main$force_intersect
  if (is.null(force_intersect)) force_intersect <- FALSE

  #' Ensure everything has the same gene IDs; if so then we'll grab just the
  #' gene IDs from the first one to use as our keys.
  nms <- lapply(res_list, function(x) rownames(x$res))
  if (!all(sapply(nms, function(x) all(x==nms[[1]])))){
    if (!force_intersect){
      stop("Gene IDs differ between results objects. Consider using force_intersect=TRUE to only keep the common set.")
    } else {

      # Take advantage of utilities originally used for UpSet plots to get the
      # set of genes found in all lists
      warning("Keeping only the intersecting set of genes across all results lists because force_intersect=TRUE.")
      keep <- lcdbwf:::fromList.with.names(nms)
      keep <- rownames(keep[rowSums(keep) == ncol(keep),])
      res_list <- lapply(res_list, function (x) {
          x$res <- x$res[keep,]
          x}
      )
    }
  }

  keys <- rownames(res_list[[1]]$res)

  if (use_orgdb == TRUE) {
    orgdb <- lcdbwf:::get_annotation_db(config, dbtype="OrgDb")

    # Create a dataframe mapping gene IDs to the various configured columns
    lookups <- list()
    for (col in config$annotation$orgdb_columns){
      lookups[[col]] <- mapIds(orgdb, keys=keys, column=col, keytype=config$annotation$keytype, multiVal='first')
      if (col %in% config$annotation$fill){
        lookups[[col]] <- ifelse(is.na(lookups[[col]]), keys, lookups[[col]])
    }
    }
    lookups <- data.frame(lookups)

    # Use that dataframe to attach additional columns to each results object
    for (name in names(res_list)){
      res <- res_list[[name]]$res
      orig_colnames <- colnames(res)
      for (col in config$annotation$orgdb_columns){
        res[[col]] <- lookups[[col]]
      }
      res$gene <- rownames(res)
      res <- res[, c('gene', config$annotation$orgdb_columns, orig_colnames)]
      res_list[[name]]$res <- res
    }
  } else {
    # attach the genes as SYMBOLs in absence of OrgDb data
    for (name in names(res_list)){
      res <- res_list[[name]]$res
      res$gene <- rownames(res)
      res$SYMBOL <- rownames(res)
      res_list[[name]]$res <- res
    }
  }

  return(res_list)
}
