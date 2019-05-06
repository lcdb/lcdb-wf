import os, sys
sys.path.append('../../')
from textwrap import dedent
import tempfile
from snakemake.shell import shell
from lib import utils
log = snakemake.log_fmt_shell(append=True)

# Since we'll be appending the output from multiple commands to the same log,
# we want to ensure that the provided log file is empty to start
if snakemake.log:
    shell('cat /dev/null > {snakemake.log}')

java_args = snakemake.params.get('java_args', '')
samtools_merge_extra = snakemake.params.get('samtools_merge_extra', '')
markduplicates_extra = snakemake.params.get('markduplicates_extra', '')

if len(snakemake.input) == 1:
    utils.make_relative_symlink(snakemake.input[0], snakemake.output.bam)
    shell('touch {snakemake.output.metrics}')

else:

    merged = tempfile.NamedTemporaryFile(delete=False, prefix='merged', suffix='.bam').name
    merged_and_deduped = snakemake.output.bam
    if 'metrics' in snakemake.output.keys():
        metrics = snakemake.output.metrics
    else:
        metrics = snakemake.output.bam + '.metrics'

    shell('echo "tempfile created: {merged}" {log}')

    shell(
        'samtools merge '
        '-f '
        '-@ {snakemake.threads} '
        '{samtools_merge_extra} '
        '{merged} '
        '{snakemake.input} '
        '{log} '
    )
    shell(
        'picard '
        '{java_args} '
        'MarkDuplicates '
        'INPUT={merged} '
        'OUTPUT={merged_and_deduped} '
        'METRICS_FILE={metrics} '
        'REMOVE_DUPLICATES=true '
        '{markduplicates_extra} '
        '{log} '
    )

    shell('rm -v {merged} {log}')
