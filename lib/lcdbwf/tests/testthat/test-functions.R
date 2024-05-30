# Helper function to make minimal default design data
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

# Function to create design data for LRT test
make_lrt_design_data <- function() {
  lrt_design_data <- make_design_data()
  lrt_design_data$test <- 'LRT'
  lrt_design_data$reduced_design <- ~1
  return(lrt_design_data)
} # make_lrt_design_data

make_dds_list <- function() {
  # Create design data and dds object for Wald test type
  wald_design_data <- make_design_data()
  dds_wald <- make_dds(wald_design_data,
                       config=config,
                       featureCounts='featurecounts.txt',
                       parallel=config$parallel$parallel)

  lrt_design_data <- make_lrt_design_data()
  dds_lrt <- make_dds(lrt_design_data,
                      config=config,
                      featureCounts='featurecounts.txt',
                      parallel=config$parallel$parallel)

  # Create dds_list
  dds_list <- list(dds_wald=dds_wald, dds_lrt=dds_lrt)
  return(dds_list)
} # make_dds_list

