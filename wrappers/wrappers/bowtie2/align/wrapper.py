__author__ = "Ryan Dale"
__copyright__ = "Copyright 2016, Ryan Dale"
__email__ = "dalerr@niddk.nih.gov"
__license__ = "MIT"

from snakemake.shell import shell
from lcdblib.snakemake import aligners

bowtie2_extra = snakemake.params.get('bowtie2_extra', '')
samtools_view_extra = snakemake.params.get('samtools_view_extra', '')
samtools_sort_extra = snakemake.params.get('samtools_sort_extra', '')

log = snakemake.log_fmt_shell()

# Handle paired-end reads. Since snakemake automatically converts a one-element
# list to a string, here we detect single-end reads by checking if input.fastq
# is a string.
if isinstance(snakemake.input.fastq, str):
    fastqs = '-U {0} '.format(snakemake.input.fastq)
else:
    assert len(snakemake.input.fastq) == 2
    fastqs = '-1 {0} -2 {1} '.format(*snakemake.input.fastq)

prefix = aligners.prefix_from_bowtie2_index(snakemake.input.index)

output_prefix = snakemake.output.bam.replace('.bam', '')

shell(
    "bowtie2 "
    "-x {prefix} "
    "{fastqs} "
    "--threads {snakemake.threads} "
    "{bowtie2_extra} "
    "-S {output_prefix}.sam "
    "{log}"
)

# bowtie2 outputs SAM format so we convert to BAM here.
shell(
    "samtools view -Sb "
    "{samtools_view_extra} "
    "{output_prefix}.sam "
    "> {output_prefix}.tmp.bam && rm {output_prefix}.sam"
)

# sort the BAM and clean up
shell(
    "samtools sort "
    "-o {snakemake.output.bam} "
    "{samtools_sort_extra} "
    "-O BAM "
    "{output_prefix}.tmp.bam "
    "&& rm {output_prefix}.tmp.bam "
)
