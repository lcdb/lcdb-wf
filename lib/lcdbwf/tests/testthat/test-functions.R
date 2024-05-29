# Helper function to make minimal default design data. design_data is an argument and
# object of type list that is passed to make_dds()
make_design_data <- function() {
  lst <- list(
  # Create the sample table
  sampletable = read.table("../../../../workflows/rnaseq/config/sampletable.tsv",
                           sep="\t",
                           header=TRUE),
  design = ~ group
  ) # lst
  lst$sampletable$group <- as.factor(lst$sampletable$group)
  return(lst)
} # make_default_wald_design_data

