#!/usr/bin/env python

import os
from tempfile import NamedTemporaryFile, gettempdir
from snakemake.shell import shell

# All wrappers must be able to handle an optional params.extra.
extra = snakemake.params.get('extra', '')
# This lets us handle whether to write to a log file or to write to stdout.
# See snakemake.script.log_fmt_shell for details.
if snakemake.log:
    snakemake.log = os.path.realpath(str(snakemake.log))

log = snakemake.log_fmt_shell(stdout=False)

# Get directories that I need to move between
cwd = os.getcwd()
tmpdir = gettempdir()

# tin uses the name of the BAM to create outputs. In order to write outputs to
# tmp I need to copy the BAM over to tmp.
bam = NamedTemporaryFile(suffix='.bam').name
bed = NamedTemporaryFile(suffix='.bed').name
name = bam.rstrip('.bam')

shell(
    'cp {snakemake.input.bam} {bam} '
    '&& cp {snakemake.input.bam}.bai {bam}.bai '
    '&& cp {snakemake.input.bed} {bed}')

os.chdir(tmpdir)
shell(
    'tin.py '
    '-i {bam} '
    '-r {bed} '
    '{extra} '
    '{log}')

# Cleanup 1
shell(
    'rm {bam} '
    '&& rm {bam}.bai '
    '&& rm {bed}')

# Move outputs
os.chdir(cwd)
shell(
    'cp {name}.tin.xls {snakemake.output.table} '
    '&& cp {name}.summary.txt {snakemake.output.summary} '
    '&& rm {name}.tin.xls '
    '&& rm {name}.summary.txt')
