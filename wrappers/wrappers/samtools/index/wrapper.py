from snakemake.shell import shell
extra = snakemake.params.get('extra', '')
log = snakemake.log_fmt_shell()
shell("samtools index {extra} {snakemake.input.bam} {snakemake.output.bai}")
