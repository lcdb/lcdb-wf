__author__ = "Ryan Dale"
__copyright__ = "Copyright 2016, Ryan Dale"
__email__ = "dalerr@niddk.nih.gov"
__license__ = "MIT"

from snakemake.shell import shell

extra = snakemake.params.get('extra', '')
log = snakemake.log_fmt_shell()

shell(
    "featureCounts "
    "-T {snakemake.threads} "
    "{extra} "
    "-a {snakemake.input.annotation} "
    "-o {snakemake.output.counts} "
    "{snakemake.input.bam} "
    "{log} ")
