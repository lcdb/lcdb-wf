# ------------------------------------------------------------------------------
# Functions for working with annotations and OrgDbs
# ------------------------------------------------------------------------------

get_annotation_hub <- function(config, localHub=NULL, force=NULL, cache=NULL){

  if (missing(localHub)) localHub <- config$annotation$localHub
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

  ah <- AnnotationHub(
    hub=getAnnotationHubOption('URL'),
    proxy=proxy,
    localHub=localHub,
    cache=cache
  )
  return(ah)
}


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
get_annotation_db <- function(config, dbtype, genus_species=NULL, orgdb_key_override=NULL, txdb_key_override=NULL, cache=NULL){

    if (missing(genus_species)) genus_species <- config$annotation$genus_species
    if (missing(cache)) cache <- config$annotation$hub_cache
    if (missing(orgdb_key_override)) orgdb_key_override <- config$annotation$orgdb_key_override
    if (missing(txdb_key_override)) txdb_key_override <- config$annotation$txdb_key_override

    ah <- get_annotation_hub(config, cache=cache)

    # If an override was provided, immediately return the corresponding db.
    if (!is.null(orgdb_key_override) & dbtype == "OrgDb") return(ah[[orgdb_key_override]])
    if (!is.null(txdb_key_override) & dbtype == "TxDb") return(ah[[txdb_key_override]])

    ah.query <- query(ah, c(genus_species, dbtype))

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
      arrange(desc(rdatadateadded)) %>%
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
#' (e.g., fill in missing symbols with Ensembl IDs)
#'
#' @param res_list List of DESeqResults objects
#' @param config Full config object, at least containing config$orgdb
#'
#' @return List of same results objects, but each one with additional columns
#'   attached as specified in the config
attach_extra <- function(res_list, config, force_intersect){

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
      keep <- fromList.with.names(nms)
      keep <- rownames(keep[rowSums(keep) == ncol(keep),])
      res_list <- lapply(res_list, function (x) {
          x$res <- x$res[keep,]
          x}
      )
    }
  }

  keys <- rownames(res_list[[1]]$res)

  orgdb <- lcdbwf::get_annotation_db(config, dbtype="OrgDb")

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
  return(res_list)
}
