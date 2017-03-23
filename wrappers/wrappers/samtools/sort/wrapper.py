import os
from snakemake.shell import shell

prefix = os.path.splitext(snakemake.output[0])[0]
log = snakemake.log_fmt_shell()
extra = snakemake.params.get('extra', '')

shell(
    "samtools sort "
    "{extra} "
    "-@ {snakemake.threads} "
    "-o {snakemake.output.bam} "
    "-T {prefix} "
    "{snakemake.input.bam} "
)
