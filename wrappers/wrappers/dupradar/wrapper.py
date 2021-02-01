import tempfile
from snakemake.shell import shell
import os, sys
sys.path.append(os.path.abspath('../..'))
from lib import helpers

extra = snakemake.params.get('extra', '')
try:
    log = snakemake.log
except AttributeError:
    log = None

stranded = snakemake.params.get('stranded', False)
try:
    stranded_int = {False: 0, True: 1, 'reverse': 2}[stranded]
except KeyError:
    raise ValueError('"stranded" must be True|False|"reverse"')

paired = snakemake.params.get('paired', False)
try:
    paired_bool= {True: 'TRUE', False: 'FALSE'}[paired]
except KeyError:
    raise ValueError('"paired" must be True or False')

tempdir = tempfile.mkdtemp()

# To avoid issues with png() related to X11 and cairo, we can use bitmap() instead.
# (thanks
# http://stackoverflow.com/questions/24999983/
# r-unable-to-start-device-png-capabilities-has-true-for-png
# #comment52353278_25064603 )

script = """
library(dupRadar)
bam <- "{snakemake.input.bam}"
gtf <- "{snakemake.input.annotation}"
dm <- analyzeDuprates(bam, gtf, {stranded_int}, {paired_bool}, {snakemake.threads}, tmpDir = "{tempdir}")

dm$mhRate <- (dm$allCountsMulti - dm$allCounts) / dm$allCountsMulti
bitmap(file="{snakemake.output.multimapping_histogram}")
hist(dm$mhRate, breaks=50, main=basename(bam),
    xlab="Multimapping rate per gene", ylab="Frequency")
dev.off()

bitmap(file="{snakemake.output.density_scatter}")
duprateExpDensPlot(dm, main=basename(bam))
dev.off()

bitmap(file="{snakemake.output.expression_histogram}")
expressionHist(dm)
dev.off()

bitmap(file="{snakemake.output.expression_boxplot}")
par(mar=c(10,4,4,2)+.1)
duprateExpBoxplot(dm, main=basename(bam))
dev.off()

bitmap(file="{snakemake.output.expression_barplot}")
readcountExpBoxplot(dm)
dev.off()

write.table(dm, file="{snakemake.output.dataframe}", sep="\\t")

# The following is from
# https://github.com/ewels/NGI-RNAseq/blob/master/bin/dupRadar.r

fit <- duprateExpFit(DupMat=dm)
df <- data.frame(intercept=as.numeric(fit$intercept), slope=c(fit$slope))
cat("# dupRadar model params\\n", file="{snakemake.output.model}")
write.table(df, file="{snakemake.output.model}", sep="\\t", append=TRUE, row.names=FALSE)

# Get numbers from dupRadar GLM
curve_x <- sort(log10(dm$RPK))
curve_y = 100*predict(fit$glm, data.frame(x=curve_x), type="response")
# Remove all of the infinite values
infs = which(curve_x %in% c(-Inf,Inf))
curve_x = curve_x[-infs]
curve_y = curve_y[-infs]
# Reduce number of data points
curve_x <- curve_x[seq(1, length(curve_x), 10)]
curve_y <- curve_y[seq(1, length(curve_y), 10)]
# Convert x values back to real counts
curve_x = 10^curve_x
# Write to file
write.table(
  cbind(curve_x, curve_y),
  file="{snakemake.output.curve}",
  quote=FALSE, row.names=FALSE
)
""".format(**locals())

tmp = tempfile.NamedTemporaryFile(delete=False).name
helpers.rscript(script, tmp, log=log)
shell("rm -r {tempdir}")
