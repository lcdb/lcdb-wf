import os
import math
import tempfile
from snakemake import shell

log = snakemake.log_fmt_shell()
logfile = None
#extra = snakemake.params.get('extra', '')

outdir, basebed = os.path.split(snakemake.output.bed)
label = snakemake.params.block['label']
#extra = snakemake.params.block.get('extra', '')

ip_bdg = snakemake.input.get("ip_bdg")
lambda_bdg = snakemake.input.get("lambda_bdg")
d_estimate = snakemake.params.get("d_estimate")
read_length = snakemake.params.get("read_length")
bed = snakemake.output.get("bed")

qvalue = snakemake.params.block.get('qvalue')
if qvalue is None:
    qvalue = 0.05
else:
    qvalue = float(qvalue[0])

qvalue = -1 * math.log(qvalue, 10)

if ip_bdg is None:
    raise ValueError("macs2/bdgpeakcall requires input.ip_bdg")
if lambda_bdg is None:
    raise ValueError("macs2/bdgpeakcall requires input.lambda_bdg")
if d_estimate is None:
    raise ValueError("macs2/bdgpeakcall requires input.d_estimate")
if read_length is None:
    raise ValueError("macs2/bdgpeakcall requires input.read_length")
if bed is None:
    raise ValueError("macs2/bdgpeakcall requires output.bed")

with open(d_estimate) as f:
    d_estimate = f.readline().split()[2].strip()

with open(read_length) as f:
    read_length = f.readline().strip()

qvalue_intermediate = tempfile.NamedTemporaryFile().name

cmds = (
    'macs2 bdgcmp '
    '-t {ip_bdg} '
    '-c {lambda_bdg} '
    '-m qpois -o {qvalue_intermediate} '
)
shell(cmds + '{log}')

cmds = (
    'macs2 bdgpeakcall '
    '-i {qvalue_intermediate} '
    '-c {qvalue} '
    '-l {d_estimate} '
    '-g {read_length} '
    '-o {bed} '
)
shell(cmds + '{log}')

bed_cleaning_intermediate = tempfile.NamedTemporaryFile().name
# Fix the output file so that it doesn't have negative numbers and so it fits
# inside the genome
shell(
    """awk -F "\\t" '{{OFS="\\t"; print $1, "0", $2}}' """
    "{snakemake.input.chromsizes} "
    "> {bed_cleaning_intermediate}"
)
unsorted_bed_intermediate = tempfile.NamedTemporaryFile().name
shell(
    "export LC_COLLATE=C; "
    """awk -F "\\t" '{{OFS="\\t"; if (($2>0) && ($3>0)) print $0}}' {bed} | """
    "bedtools intersect -a - -b {bed_cleaning_intermediate} > {unsorted_bed_intermediate} "
    "&& bedSort {unsorted_bed_intermediate} {bed}"
)
shell("rm {0} {1} {2}".format(qvalue_intermediate,
                              bed_cleaning_intermediate,
                              unsorted_bed_intermediate))
