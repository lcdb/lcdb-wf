#!/usr/bin/env python

from snakemake.shell import shell

# All wrappers must be able to handle an optional params.extra.
extra = snakemake.params.get('extra', '')


# This lets us handle whether to write to a log file or to write to stdout.
# See snakemake.script.log_fmt_shell for details.
log = snakemake.log_fmt_shell()


# This demo shows how to handle paired-end and single-end input data as two
# different cases, depending on whether the rule's input included an "R2" key
# or not.
paired_end = (
    'R1' in snakemake.input.keys() and
    'R2' in snakemake.input.keys()
)

if paired_end: 
    shell('cp {snakemake.input.R1} {snakemake.output.R1}')
    shell('cp {snakemake.input.R2} {snakemake.output.R2}')

else:
    shell("cp {snakemake.input} {snakemake.output} {log}")
