from snakemake.shell import shell
log = snakemake.log_fmt_shell()
extra = snakemake.params.get('extra', '')
java_args = snakemake.params.get('java_args', '')
shell(
    'picard '
    '{java_args} '
    'CollectRnaSeqMetrics '
    'REF_FLAT={snakemake.input.refflat} '
    'INPUT={snakemake.input.bam} '
    'OUTPUT={snakemake.output.metrics} '
    '{extra} '
    '{log} '
)
