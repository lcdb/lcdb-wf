__author__ = "Ryan Dale"
__copyright__ = "Copyright 2016, Ryan Dale"
__email__ = "dalerr@niddk.nih.gov"
__license__ = "MIT"

from tempfile import NamedTemporaryFile
from snakemake.shell import shell
from lcdblib.snakemake import aligners

hisat2_extra = snakemake.params.get('hisat2_extra', '')
samtools_view_extra = snakemake.params.get('samtools_view_extra', '')
samtools_sort_extra = snakemake.params.get('samtools_sort_extra', '')

log = snakemake.log_fmt_shell()

# Hisat2 can either align using provided FASTQ files or given an SRA identifier.
if snakemake.input.get('fastq', ''):
    # Handle paired-end reads. Since snakemake automatically converts a one-element
    # list to a string, here we detect single-end reads by checking if input.fastq
    # is a string.
    if isinstance(snakemake.input.fastq, str):
        fastqs = '-U {0} '.format(snakemake.input.fastq)
    elif isinstance(snakemake.input.fastq, list) & (len(snakemake.input.fastq) == 1):
        fastqs = '-U {0} '.format(snakemake.input.fastq[0])
    elif isinstance(snakemake.input.fastq, list) & (len(snakemake.input.fastq) == 2):
        fastqs = '-1 {0} -2 {1} '.format(*snakemake.input.fastq)
    else:
        raise ValueError("You must provide a list of string to input.fastq")
else:
    fastqs = ''
    try:
        assert '--sra-acc' in hisat2_extra
    except:
        raise ValueError("You must provide a FASTQ file or use the --sra-acc option.")

prefix = aligners.prefix_from_hisat2_index(snakemake.input.index)

# Create temporary files to store intermediates. Will use $TMDPIR if set.
sam = NamedTemporaryFile(suffix='.sam', delete=False).name
bam = NamedTemporaryFile(suffix='.bam', delete=False).name
sort_bam = NamedTemporaryFile(suffix='.sort.bam', delete=False).name

shell(
    "hisat2 "
    "-x {prefix} "
    "{fastqs} "
    "--threads {snakemake.threads} "
    "{hisat2_extra} "
    "-S {sam} "
    "{log}"
)

# hisat2 outputs SAM format so we convert to BAM here.
shell(
    "samtools view -Sb "
    "{samtools_view_extra} "
    "{sam} "
    "> {bam} && rm {sam}"
)

# sort the BAM and clean up
shell(
    "samtools sort "
    "-o {sort_bam} "
    "{samtools_sort_extra} "
    "-O BAM "
    "{bam} "
    "&& cp {sort_bam} {snakemake.output.bam} "
    "&& rm {bam} "
    "&& rm {sort_bam} "
)
