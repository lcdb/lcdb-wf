#!/usr/bin/env python

import os
from tempfile import TemporaryDirectory
from snakemake.shell import shell

# All wrappers must be able to handle an optional params.extra.
extra = snakemake.params.get('extra', '')

# This lets us handle whether to write to a log file or to write to stdout.
# See snakemake.script.log_fmt_shell for details.
if snakemake.log:
    snakemake.log = os.path.realpath(str(snakemake.log))

log = snakemake.log_fmt_shell()

# Get directories that I need to move between
cwd = os.getcwd()

# Get filenames
bam = os.path.basename(snakemake.input.bam)
bai = os.path.basename(snakemake.input.bai)
bed = os.path.basename(snakemake.input.bed)

with TemporaryDirectory() as tmpdir:
    # Copy BAMS to tmpdir
    shell(
            'cp {snakemake.input.bam} {tmpdir}/{bam} '
            '&& cp {snakemake.input.bai} {tmpdir}/{bai}'
            '&& cp {snakemake.input.bed} {tmpdir}/{bed} '
            )

    # Move to tmpdir
    os.chdir(tmpdir)

    # Run program
    shell(
        'geneBody_coverage.py '
        '-i {bam} '
        '-o tmp '
        '-r {bed} '
        '{extra} '
        '{log}'
        )

    # Move to tmpdir
    os.chdir(cwd)
    shell(
        'cp {tmpdir}/tmp.geneBodyCoverage.r {snakemake.output.r} '
        '&& cp {tmpdir}/tmp.geneBodyCoverage.txt {snakemake.output.txt} '
        '&& cp {tmpdir}/tmp.geneBodyCoverage.curves.* {snakemake.output.img}'
        )
