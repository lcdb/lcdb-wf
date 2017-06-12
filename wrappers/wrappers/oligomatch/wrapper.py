#!/usr/bin/env python

from snakemake.shell import shell
log = snakemake.log_fmt_shell()
shell(
    "oligoMatch "
    "{snakemake.input.oligos} "
    "{snakemake.input.sequence} "
    "{snakemake.output} "
    "{log}"
)
