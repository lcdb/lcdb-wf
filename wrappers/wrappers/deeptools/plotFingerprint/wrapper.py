import os
from snakemake.shell import shell
from tempfile import NamedTemporaryFile

# Pull in parameters
extra = snakemake.params.get('extra', '')

# Make log
log = snakemake.log_fmt_shell()

if '--outQualityMetrics' not in extra and 'metrics' in snakemake.output.keys():
    extra += ' --outQualityMetrics {snakemake.output.metrics} '.format(**locals())

if '--plotFile' not in extra and 'plot' in snakemake.output.keys():
    extra += ' --plotFile {snakemake.output.plot} '.format(**locals())

if '--outRawCounts' not in extra and 'raw_counts' in snakemake.output.keys():
    extra += ' --outRawCounts {snakemake.output.raw_counts} '.format(**locals())

control = ""
if 'control' in snakemake.input.keys():
    control = ' --JSDsample {snakemake.input.control} '.format(**locals())

shell(
    'plotFingerprint '
    '--bamfiles {snakemake.input.bams} '
    '{control} '
    '-p {snakemake.threads} '
    '{extra} '
    '{log} '
)
