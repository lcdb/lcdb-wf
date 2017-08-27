import os
from snakemake.shell import shell
from lcdblib.utils import utils

prefix = os.path.splitext(snakemake.output[0])[0]
log = snakemake.log_fmt_shell()
extra = snakemake.params.get('extra', '')
if len(snakemake.input) == 0:
    raise ValueError("No inputs provided")
if len(snakemake.input) == 1:
    utils.make_relative_symlink(snakemake.input[0], snakemake.output[0])
else:
    shell(
        "samtools merge "
        "{extra} "
        "-@ {snakemake.threads} "
        "{snakemake.output} "
        "{snakemake.input} "
    )
