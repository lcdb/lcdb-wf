import os
import glob
from snakemake import shell

log = snakemake.log_fmt_shell()
logfile = None
extra = snakemake.params.get('extra', '')

outdir, basebed = os.path.split(snakemake.output.bed)
label = snakemake.params.block['label']
extra = snakemake.params.block.get('extra', '')

# `-c` has to be skipped if no control is provided
# if os.path.isfile(snakemake.input.control):
if len(snakemake.input.control) > 0:
    arguments = '-c {snakemake.input.control} '
else:
    arguments = ''
# Add `--guess-bampe` if input dataset is paired-end
if snakemake.params.is_paired:
    arguments += '--guess-bampe '

# cmds = (
#     'epic2 ' + arguments +
#     '-t {snakemake.input.ip} '
#     '--chromsizes {snakemake.input.chromsizes} | '
#     'sort -k1,1 -k2,2n > {label}.tmp.bed ' + extra
# )

shell(
    'epic2 ' + arguments + extra +
    '-t {snakemake.input.ip} '
    '--chromsizes {snakemake.input.chromsizes} | '
    'sort -k1,1 -k2,2n > {label}.tmp.bed '
)

# Fix the output file so that it doesn't have negative numbers and so it fits
# inside the genome
shell(
    """awk -F "\\t" '{{OFS="\\t"; print $1, "0", $2}}' """
    "{snakemake.input.chromsizes} "
    "> {label}.tmp.genome"
)
shell(
    "export LC_COLLATE=C; "
    """awk -F "\\t" '{{OFS="\\t"; if (($2>0) && ($3>0)) print $0}}' {label}.tmp.bed | """
    "bedtools intersect -a - -b {label}.tmp.genome > {snakemake.output.bed}.tmp "
    "&& bedSort {snakemake.output.bed}.tmp {snakemake.output.bed} "
    "&& rm {label}.tmp.bed {snakemake.output.bed}.tmp {label}.tmp.genome"
)
