from textwrap import dedent
import tempfile
from snakemake.shell import shell
log = snakemake.log_fmt_shell(append=True)

# Since we'll be appending the output from multiple commands to the same log,
# we want to ensure that the provided log file is empty to start
if snakemake.log:
    shell('cat /dev/null > {snakemake.log}')

java_args = snakemake.params.get('java_args', '')
keep_tempfiles = snakemake.params.get('keep_tempfiles', False)

registered_for_deletion = [
    snakemake.output.bed + '.tmp',
    snakemake.output.bed + '.tmp.genome',
]


def merge_and_dedup(bams):
    """
    spp only handles one replicate at a time. To support pooled samples, we
    merge and remove duplicates, storing the result in a tempfile.

    If only one item is provided, return it immediately
    """

    if len(bams) == 1:
        return bams

    merged = tempfile.NamedTemporaryFile(delete=False, prefix='merged', suffix='.bam').name
    merged_and_deduped = tempfile.NamedTemporaryFile(delete=False, prefix='merged_and_duped', suffix='.bam').name
    metrics = tempfile.NamedTemporaryFile(delete=False, prefix='metrics', suffix='.txt').name

    shell('echo "tempfiles created by merge_and_dedup: {merged} {merged_and_deduped} {metrics}" {log}')

    if not keep_tempfiles:
        registered_for_deletion.extend([merged, merged_and_deduped, metrics])

    bams = ' '.join(bams)
    shell(
        'samtools merge '
        '-f '
        '-@ {snakemake.threads} '
        '{merged} '
        '{bams} '
        '{log} '
    )
    shell(
        'picard '
        '{java_args} '
        'MarkDuplicates '
        'INPUT={merged} '
        'OUTPUT={merged_and_deduped} '
        'METRICS_FILE={metrics} '
        'REMOVE_DUPLICATES=true '
        '{log} '
    )
    return merged_and_deduped


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

DEFAULTS = {
    # srange controls the range of lags over which to calculate cross-correlation
    'srange': (50, 500),
    # bins controls how the binding characteristics will be binned
    'bins': 5,
    # enable/disable the remove.tag.anomalies step
    'remove_anomalies': False,
    # false discovery rate when calling peaks
    'fdr': 0.05,
    # window half-size. Used if binding.characteristics is NA.
    'whs': 500,
    # Z threshold used when adding broad regions.
    'zthr': 3,
    # bandwith for smoothing WIG file
    'bandwidth': 200,
    # step for smoothing WIG file
    'step': 100,
    # Set to False to disable the filtering of large regions with high input signal
    'tecfilter': True,
}

params = {}
for k, v in DEFAULTS.items():
    v = extra.get(k, v)
    if isinstance(v, bool):
        v = Rbool(v)
    params[k] = v

# ----------------------------------------------------------------------------

# R_template is incrementally built up so that we can intersperse comments and
# to keep things better organized. It will be filled in with `**locals()` at
# the end.

ip = merge_and_dedup(snakemake.input.ip)
control = merge_and_dedup(snakemake.input.control)


R_template = """
library(spp)
chip.data <- read.bam.tags("{ip}")
input.data <- read.bam.tags("{control}")
"""


#
R_template += """
for (chrom in names(chip.data$tags)){{
    if (length(chip.data$tags[[chrom]]) < 10){{
        print(paste("Chromosome", chrom, "has <10 reads; removing from analysis"))
        chip.data$tags[[chrom]] <- NULL
        chip.data$quality[[chrom]] <- NULL
        input.data$tags[[chrom]] <- NULL
        input.data$quality[[chrom]] <- NULL
    }}
}}
"""

# Use configured srange and bins, if provided. `accept.all.tags=TRUE` is
# hard-coded since we were getting errors if FALSE.
R_template += """
binding.characteristics <- get.binding.characteristics(
  chip.data,
  srange=c({params[srange][0]}, {params[srange][1]}),
  bin={params[bins]},
  accept.all.tags=TRUE,
  remove.tag.anomalies={params[remove_anomalies]}
)
"""

R_template += """
# Extract info from binding characteristics
tag.shift <- round(binding.characteristics$peak$x/2)
detection.window.halfsize <- binding.characteristics$whs
if (!is.finite(detection.window.halfsize)){{
  detection.window.halfsize <- {params[whs]}
}}
"""

R_template += """
# Reset data to tags, and remove any chromosomes with no data.
# (tags is a list, names are chromosomes and values are integer vectors)

chip.data <- chip.data$tags
input.data <- input.data$tags

chip.data[sapply(chip.data, is.null)] <- NULL
input.data[sapply(input.data, is.null)] <- NULL
"""


if 'smoothed_enrichment_mle' in snakemake.output.keys():
    R_template += dedent("""
    smoothed.enrichment.estimate <- get.smoothed.enrichment.mle(
      chip.data,
      input.data,
      bandwidth={params[bandwidth]},
      step={params[step]},
      tag.shift=tag.shift)
    writewig(
      smoothed.enrichment.estimate,
      "{snakemake.output.smoothed_enrichment_mle}",
      feature=""
    )
    """)

if 'enrichment_estimates' in snakemake.output.keys():
    R_template += dedent("""
    enrichment.estimates <- get.conservative.fold.enrichment.profile(
        chip.data, input.data, fws=500, step=100, alpha=0.01
    )
    writewig(enrichment.estimates, "{snakemake.output.enrichment_estimates}", feature="")
    rm(enrichment.estimates)
    """)

R_template += """
# Get peaks
bp <- find.binding.positions(
  signal.data=chip.data,
  control.data=input.data,
  fdr={params[fdr]},
  whs=detection.window.halfsize,
  tec.filter={params[tecfilter]}
)
"""

R_template += """
# Add broad regions to peaks
bp <- add.broad.peak.regions(
  chip.data,
  input.data,
  bp,
  window.size=detection.window.halfsize,
  z.thr={params[zthr]}
)
write.narrowpeak.binding(bp, "{snakemake.output.bed}.tmp")
"""

# Save image for later introspection or debugging
if 'rdata' in snakemake.output.keys():
    R_template += dedent("""
    save.image("{snakemake.output.rdata}")
    """)

# write the filled-in template to the output directory for later debugging
script_filename = snakemake.output.bed + '.R'
with open(script_filename, 'w') as fout:
    fout.write(R_template.format(**locals()))

# Run it
shell('Rscript {script_filename} {log}')

# Fix the output file so that it doesn't have negative numbers and so it fits
# inside the genome
shell(
    """awk -F "\\t" '{{OFS="\\t"; print $1, "0", $2}}' """
    "{snakemake.input.chromsizes} "
    "> {snakemake.output.bed}.tmp.genome"
)
shell(
    "sort -k1,1 -k2,2n {snakemake.output.bed}.tmp | "
    """awk -F "\\t" '{{OFS="\\t"; if (($2>0) && ($3>0)) print $0}}' | """
    "bedtools intersect -a - -b {snakemake.output.bed}.tmp.genome > {snakemake.output.bed}"
)

# SPP's writewig() adds a header and is space-separated, so this turns it into
# a proper bedGraph file ready for conversion to bigwig.
if 'enrichment_estimates' in snakemake.output.keys():
    shell('grep -v "track" {snakemake.output.enrichment_estimates} '
          '| sed "s/ /\\t/g" > {snakemake.output.enrichment_estimates}.tmp '
          '&& mv {snakemake.output.enrichment_estimates}.tmp '
          '{snakemake.output.enrichment_estimates}')

if 'smoothed_enrichment_mle' in snakemake.output.keys():
    shell('grep -v "track" {snakemake.output.smoothed_enrichment_mle} '
          '| sed "s/ /\\t/g" > {snakemake.output.smoothed_enrichment_mle}.tmp '
          '&& mv {snakemake.output.smoothed_enrichment_mle}.tmp '
          '{snakemake.output.smoothed_enrichment_mle}')

for fn in registered_for_deletion:
    shell('rm -v {fn} {log}')
