import os
from snakemake.shell import shell
log = snakemake.log_fmt_shell()


def Rbool(x):
    """
    Convert to R boolean string used to fill in a template
    """
    if x:
        return 'TRUE'
    return 'FALSE'

# ----------------------------------------------------------------------------
# DEFAULTS
#
extra = snakemake.params.block.get('extra', {})

# srange controls the range of lags over which to calculate cross-correlation
srange = extra.get('srange', (50, 500))

# bins controls how the binding characteristics will be binned
bins = extra.get('bins', 5)

# enable/disable the remove.tag.anomalies step
remove_anomalies = Rbool(extra.get('remove_anomalies', True))

# false discovery rate when calling peaks
fdr = extra.get('fdr', 0.05)

# window half-size. Used if binding.characteristics is NA.
whs = extra.get('whs', 500)

# Z threshold used when adding broad regions.
zthr = extra.get('zthr', 3)

# bandwith for smoothing WIG file
bandwidth = extra.get('bandwidth', 200)

# step for smoothing WIG file
step = extra.get('step', 100)

# ----------------------------------------------------------------------------

# R_template is incrementally built up so that we can intersperse comments and
# to keep things better organized. It will be filled in with `**locals()` at
# the end.
R_template = """
library(spp)
chip.data <- read.bam.tags("{snakemake.input.ip}")
input.data <- read.bam.tags("{snakemake.input.control}")
"""

# Use configured srange and bins, if provided. `accept.all.tags=TRUE` is
# hard-coded since we were getting errors if FALSE.
R_template += """
binding.characteristics <- get.binding.characteristics(
  chip.data,
  srange=c({srange[0]}, {srange[1]}),
  bin={bins},
  accept.all.tags=TRUE,
  remove.tag.anomalies={remove_anomalies}
)
"""

if 'smoothed_enrichment' in snakemake.output:
    R_template += """
    smoothed.enrichment.estimate <- get.smoothed.enrichment.mle(
      chip.data,
      input.data,
      bandwidth={bandwidth},
      step={step},
      tag.shift=tag.shift)
    writewig(
      smoothed.enrichment.estimate,
      "{snakemake.output.smoothed_enrichment}"
    )
    """

R_template += """
# Extract info from binding characteristics
tag.shift <- round(binding.characteristics$peak$x/2)
detection.window.halfsize <- binding.characteristics$whs
if (!is.finite(detection.window.halfsize)){{
  detection.window.halfsize <- {whs}
}}
"""

R_template += """
# Reset data to tags, and remove any chromosomes with no data

chip.data <- chip.data$tags
input.data <- input.data$tags

chip.data[sapply(chip.data, is.null)] <- NULL
input.data[sapply(input.data, is.null)] <- NULL
"""

R_template += """
# Get peaks
bp <- find.binding.positions(
  signal.data=chip.data,
  control.data=input.data,
  fdr={fdr},
  whs=detection.window.halfsize
)
"""

# peaks
R_template += """
write.narrowpeak.binding(bp, "{snakemake.output.bed}")
bp <- add.broad.peak.regions(
  chip.data,
  input.data,
  bp,
  window.size=detection.window.halfsize,
  z.thr={zthr}
)
write.narrowpeak.binding(bp, paste0("{snakemake.output.bed}", ".broadPeak"))
"""

# Save image for later introspection or debugging
image_filename = os.path.join(os.path.dirname(snakemake.output.bed), 'image.RData')
R_template += """
save.image("{image_filename}")
"""

# write the filled-in template to the output directory for later debugging
script_filename = os.path.join(os.path.dirname(snakemake.output.bed), 'script.R')
with open(script_filename, 'w') as fout:
    fout.write(R_template.format(**locals()))

# Run it
shell('Rscript {script_filename} {log}')
