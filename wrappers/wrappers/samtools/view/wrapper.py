import os
from snakemake.shell import shell

prefix = os.path.splitext(snakemake.output[0])[0]
log = snakemake.log_fmt_shell()
extra = snakemake.params.get('extra', '')

shell(
    "samtools view "
    "{extra} "
    "{snakemake.input} "
    "> {snakemake.output} "
    "{log} "
)
