import os
from snakemake.shell import shell
from tempfile import NamedTemporaryFile

# Pull in parameters
extra = snakemake.params.get('extra', '')

# Make log
log = snakemake.log_fmt_shell()

# Copy files over to TMPDIR
tmp_in = NamedTemporaryFile(suffix='.bam').name
shell(
    'cp {snakemake.input.bam} {tmp_in}'
    '&& cp {snakemake.input.bai} {tmp_in}.bai'
)

# Run bamCoverage
tmp_out = NamedTemporaryFile(suffix=os.path.splitext(snakemake.output[0])[1]).name
shell(
    'bamCoverage '
    '--bam {tmp_in} '
    '-o {tmp_out} '
    '-p {snakemake.threads} '
    '{extra} '
    '{log} '
    '&& mv {tmp_out} {snakemake.output[0]} '
    '&& rm {tmp_in}'
)
