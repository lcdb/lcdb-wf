import os
from snakemake.shell import shell

extra = snakemake.params.get('extra', '')
log = snakemake.log_fmt_shell()

outdir = [os.path.dirname(i) for i in snakemake.output]

if len(set(outdir)) != 1:
    raise ValueError("Inconsistent output directories: {0}".format(set(outdir)))

outdir = outdir[0]
print(os.path.abspath(outdir))
cmd = (
    "kallisto quant "
    "--index {snakemake.input.index} "
    "--output-dir {outdir} "
    "--threads {snakemake.threads} "
    "{extra} "
    "{snakemake.input.fastq} "
    "{log} ".format(**locals()))
print(cmd)
shell(cmd)
