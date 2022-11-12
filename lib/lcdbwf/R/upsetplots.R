#-------------------------------------------------------------------------------
# Functions for building UpSet plots out of collections of DESeqResults objects
#-------------------------------------------------------------------------------


#' Plot UpSet plots
#'
#' @param res_list "reslist" object
#' @param label_column When creating output files, use this column to label genes
plot_upsets <- function(res_list, label_column=NULL, alpha=0.1, lfc_thresh=0){

  if (!is.null(label_column)){
    lookup <- data.frame(x=res_list[[1]][['res']][[label_column]])
    colnames(lookup) <- c(label_column)
    rownames(lookup) <- rownames(res_list[[1]][['res']])
  }

  sel.list <- lapply(
    res_list,
    function(x) {
      list(
        up=lcdbwf:::get_sig(x$res, alpha=alpha, lfc_thresh=lfc_thresh, "up", return_type="rownames"),
        down=lcdbwf:::get_sig(x$res, alpha=alpha, lfc_thresh=lfc_thresh, "down", return_type="rownames"),
        changed=lcdbwf:::get_sig(x$res, alpha=alpha, lfc_thresh=lfc_thresh, "changed", return_type="rownames")
      )
    }
  )

  # Transpose the list so that instead of contrast -> direction it's
  # direction -> contrast
  sel.list.transformed <- purrr::transpose(sel.list)

  for (sel.name in names(sel.list.transformed)){
      lcdbwf:::mdcat("## UpSet plot: ", sel.name)

      # Only consider contrasts with sig genes in this selection
      ll <- purrr::keep(
          sel.list.transformed[[sel.name]],
          function(x) length(x) > 0
      )

      if (length(ll) <= 1){
          lcdbwf:::mdcat('not enough contrasts with upregulated genes')
          next
      }

      print(upset(UpSetR::fromList(ll), order.by='freq', nsets=length(ll)))
      upset.df <- lcdbwf:::fromList.with.names(ll)
      upset.df$names <- rownames(upset.df)

      first_cols <- c('names')
      if (!is.null(label_column)){
        upset.df[[label_column]] <- lookup[upset.df$names, label_column]
        first_cols <- c(first_cols, label_column)
      }

      sort_cols <- colnames(df) %>% setdiff(first_cols)
      upset.df <- upset.df %>%
          dplyr::relocate(!!!first_cols) %>%
          dplyr::arrange(!!!sort_cols)

      lcdbwf:::write.upset.plot.results(upset.df, sel.name)
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


#' Save TSVs and PDFs of UpSet plot. Saves the current plot to PDF, so should
#' be called immediately after plotting.
#'
#' @param lldf Dataframe of 1 and 0 from creating UpSet plots
#' @param label Label to use when constructing filenames
write.upset.plot.results <- function(lldf, label){
  outfile <- file.path('results', 'upset_plots', paste0(label, '.tsv'))
  dir.create(dirname(outfile), showWarnings=FALSE, recursive=TRUE)
  write.table(lldf, file=outfile, sep='\t', row.names=FALSE, quote=FALSE)
  lcdbwf:::mdcat('- [', outfile, ']', '(', outfile, ')')
  pdf.file <- file.path('results', 'upset_plots', paste0(label, '.pdf'))
  dev.copy(pdf, file=pdf.file)
  dev.off()
  lcdbwf:::mdcat('- [', pdf.file, ']', '(', pdf.file, ')')
}

