from snakemake.shell import shell
import json

extra = snakemake.params.get('extra', '')
log = snakemake.log_fmt_shell()

shell(
    'kallisto index {snakemake.input.fasta} '
    '--index {snakemake.output.index} '
    '{extra} '
    '{log} '
)
