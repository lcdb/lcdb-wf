import os
from snakemake.shell import shell
extra = snakemake.params.get('extra', '')
log = snakemake.log_fmt_shell()

outdir = [os.path.dirname(i) for i in snakemake.output]

if len(set(outdir)) != 1:
    raise ValueError("Inconsistent output directories: {0}".format(set(outdir)))

outdir = outdir[0]
if outdir == '':
    outdir = '.'

if isinstance(snakemake.input.index, str):
    index_dir = [os.path.dirname(snakemake.input.index)]
else:
    index_dir = [os.path.dirname(i) for i in snakemake.input.index]

if len(set(index_dir)) != 1:
    raise ValueError("Inconsistent index directories: {0}".format(set(index_dir)))

cmd = (
    "salmon quant "
    "--index {index_dir} "
    "--output {outdir} "
    "--threads {snakemake.threads} "
)
if 'unmatedReads' in snakemake.input.keys():
    cmd += '--unmatedReads {snakemake.input.unmatedReads} '
elif 'mates2' in snakemake.input and 'mates1' in snakemake.input.keys():
    cmd += '--mates1 {snakemake.input.mates1} --mates2 {snakemake.input.mates2} '
else:
    raise ValueError(
        "Use input key 'unmatedReads' for single-end; 'mates1' and 'mates2' for "
        "paired-end. Got: {}".format(snakemake.input.keys())
    )

cmd += (
    "{extra} "
    "{log} ".format(**locals())
)
shell(cmd)
