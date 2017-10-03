import os
from snakemake.shell import shell

outdir = os.path.dirname(snakemake.output[0])
if not outdir:
    outdir = '.'
basename = os.path.basename(snakemake.output[0])

extra = snakemake.params.get('extra', "")
log = snakemake.log_fmt_shell()

# MultiQC uses Click, which in turn complains if C.UTF-8 is not set
shell(
    'LC_ALL=en_US.UTF.8 LC_LANG=en_US.UTF-8 '
    'multiqc '
    '--quiet '
    '--outdir {outdir} '
    '--force '
    '--filename {basename} '
    '{extra} '
    '{snakemake.params.analysis_directory} '
    '{log}'
)
