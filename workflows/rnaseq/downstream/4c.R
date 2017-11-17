library(yaml)
library(argparse)
parser <- ArgumentParser()
parser$add_argument('--config', nargs=1, help='config.yaml file')
parser$add_argument('--comparison', nargs=1, help='Comparison label to run')
comparison <- 'shep-vs-gfp-brat'
config.filename <- 'config/4c-config.yml'

args <- parser$parse_args()
comparison <- args$comparison
config.filename <- args$config

config <- yaml.load_file(config.filename)
sampletable.filename <- config[['sampletable']]

sampletable <- read.table(sampletable.filename, header=TRUE, stringsAsFactors=FALSE)
dir.4c <- config[['4c_dir']]
dir.sample <- config[['sample_dir']]
config <- config[['4c']]
comparisons <- config[['comparisons']][[comparison]]
control.samples <- comparisons[['control']]
treatment.samples <- comparisons[['treatment']]
samples <- c(control.samples, treatment.samples)

# get the enzyme, bait, fraglengths, and restriction sites from the sampletable
# and make sure they are unique
meta <- list()
for (name in c('enzyme', 'bait', 'fragLen', 'fragLen2', 'restriction_site')){
    x <- unique(sampletable[sampletable$samplename %in% samples, name])
    if (length(x) != 1){
        stop(paste('More than one', name, 'specified for these samples'))
    }
    meta[[name]] <- x[1]
}
enzyme <- meta[['enzyme']]
bait <- meta[['bait']]
fragLen <- meta[['fragLen']]
fragLen2 <- meta[['fragLen2']]
restriction_site <- meta[['restriction_site']]

bait_chr <- config[['baits']][[bait]][['chrom']]
bait_coord <- config[['baits']][[bait]][['pos']]

cis_k <- config[['baits']][[bait]][['cis_k']]
nearbait_k <- config[['baits']][[bait]][['nearbait_k']]

out.base <- file.path('4cker-output', comparison)
dir.create(out.base, recursive=TRUE, showWarnings=FALSE)
enz_file <- read.table(file.path(dir.4c, paste0(enzyme, '_', fragLen, 'mer_flanking_sites_unique.bed')))
bedgraph.filename <- function (x) {
    file.path(dir.sample, x,
              paste0(enzyme, '_', fragLen, 'mer_trim', fragLen2),
              paste0(x, '_cleaned.bedgraph')
              )
}
treatment.files <- as.character(sapply(treatment.samples, bedgraph.filename))
control.files <- as.character(sapply(control.samples, bedgraph.filename))

library(R.4Cker)
both.obj <- createR4CkerObjectFromFiles(
  files=c(control.files, treatment.files),
  samples=c(control.samples, treatment.samples),
  primary_enz=restriction_site,
  bait_chr=bait_chr,
  bait_coord=bait_coord,
  bait_name=bait,
  conditions=c('control', 'treatment'),
  replicates=c(length(control.files), length(treatment.files)),
  species='dm',
  output_dir=out.base,
  enz_file=enz_file)


ensure.output <- function (obj, kind) {
    for (sample in obj@samples){
        for (suffix in c('_highinter.bed', '_lowinter.bed', '_noninter.bed')){
            fn <- file.path(
                 obj@output_dir,
                 paste0(sample, '_', kind, suffix)
            )
            if (!file.exists(fn)){file.create(fn)}
        }
    }
}


outdir <- file.path(out.base, paste0('nearbait_k', nearbait_k, '/'))
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
both.obj@output_dir <- outdir
nearbait.res <- try(nearBaitAnalysis(both.obj, k=nearbait_k))
ensure.output(both.obj, 'nearbait')
df <- nearbait.res$window_counts
names(df) <- c('chrom', 'start', 'end', 'distance', c(control.samples, treatment.samples))
pdf(file=file.path(outdir, 'splom.pdf'))
pairs(log(df[,5:ncol(df)] + 1))
dev.off()

outdir <- file.path(out.base, paste0('cis_k', cis_k, '/'))
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
both.obj@output_dir <- outdir
cis.res <- try(cisAnalysis(both.obj, k=cis_k))
ensure.output(both.obj, 'cis')
df <- cis.res$window_counts
names(df) <- c('chrom', 'start', 'end', 'distance', c(control.samples, treatment.samples))
pdf(file=file.path(outdir, 'splom.pdf'))
pairs(log(df[,5:ncol(df)] + 1))
dev.off()

outdir <- file.path(out.base, paste0('nearbait_diff_k',nearbait_k, '/'))
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
both.obj@output_dir <- outdir
differentialAnalysis(
  obj=both.obj,
  norm_counts_avg=nearbait.res$norm_counts_avg,
  windows=nearbait.res$window_counts,
  conditions=c('control', 'treatment'),
  region='nearbait',
  coordinates=NULL,
  pval=0.05)


outdir <- file.path(out.base, paste0('cis_diff_k',cis_k, '/'))
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
both.obj@output_dir <- outdir
differentialAnalysis(
  obj=both.obj,
  norm_counts_avg=cis.res$norm_counts_avg,
  windows=cis.res$window_counts,
  conditions=c('control', 'treatment'),
  region='cis',
  coordinates=NULL,
  pval=0.05)
