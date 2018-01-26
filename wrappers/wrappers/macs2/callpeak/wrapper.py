import os
import glob
from snakemake import shell

log = snakemake.log_fmt_shell()
logfile = None
extra = snakemake.params.get('extra', '')

outdir, basebed = os.path.split(snakemake.output.bed)
label = snakemake.params.block['label']
extra = snakemake.params.block.get('extra', '')

cmds = (
    'macs2 '
    'callpeak '
    '-c {snakemake.input.control} '
    '-t {snakemake.input.ip} '
    '-f BAM '
    '--outdir {outdir} '
    '--name {label} '
)
# add any per-peak-calling-run extra commands
cmds += extra

shell(cmds + ' {log}')

# Depending on whether --broad was used, we may have a "*_peaks.narrowPeak" or
# a "*_peaks.broadPeak". Figure out which it is, and fix the output to have
# non-negative values and truncate it to the genome.

narrow = glob.glob(os.path.join(outdir, '*_peaks.narrowPeak'))
broad = glob.glob(os.path.join(outdir, '*_peaks.broadPeak'))

if len(narrow) == 1:
    hit = narrow[0]
    basehit = os.path.basename(narrow[0])
elif len(broad) == 1:
    hit = broad[0]
    basehit = os.path.basename(broad[0])
else:
    raise ValueError("No narrowPeak or broadPeak found!")


# Fix the output file so that it doesn't have negative numbers and so it fits
# inside the genome
shell(
    """awk -F "\\t" '{{OFS="\\t"; print $1, "0", $2}}' """
    "{snakemake.input.chromsizes} "
    "> {hit}.tmp.genome"
)
shell(
    "sort -k1,1 -k2,2n {hit} | "
    """awk -F "\\t" '{{OFS="\\t"; if (($2>0) && ($3>0)) print $0}}' | """
    "bedtools intersect -a - -b {hit}.tmp.genome > {snakemake.output.bed}"
)
