__author__ = "Ryan Dale"
__copyright__ = "Copyright 2016, Ryan Dale"
__email__ = "dalerr@niddk.nih.gov"
__license__ = "MIT"

from snakemake.shell import shell
from lcdblib.snakemake import aligners

extra = snakemake.params.get('extra', '')
log = snakemake.log_fmt_shell()

prefix = aligners.prefix_from_hisat2_index(snakemake.output.index)
shell(
    "hisat2-build "
    "{extra} "
    "{snakemake.input.fasta} "
    "{prefix} "
    "{log} "
)
