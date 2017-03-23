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

# Copy files over to ease I/O on filesystem.
bam = NamedTemporaryFile(suffix='.bam').name
bed = NamedTemporaryFile(suffix='.bed').name
name = bam.rstrip('.bam')

shell(
    'cp {snakemake.input.bam} {bam} '
    '&& cp {snakemake.input.bed} {bed}')

os.chdir(tmpdir)
shell(
    'infer_experiment.py '
    '-i {bam} '
    '-r {bed} '
    '{extra} '
    '> {name}.txt '
    '{log}')

# Cleanup 1
shell(
    'rm {bam} '
    '&& rm {bed}')

# Move outputs
os.chdir(cwd)
shell(
    'cp {name}.txt {snakemake.output.txt} '
    '&& rm {name}.txt')
