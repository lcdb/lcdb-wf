import os
from snakemake import shell
import tempfile

log = snakemake.log_fmt_shell()
logfile = None
extra = snakemake.params.get('extra', '')

label = snakemake.params.block['label']
extra = snakemake.params.block.get('extra', '')


d_bdg = snakemake.input.d_bdg
slocal_bdg = snakemake.input.slocal_bdg
llocal_bdg = snakemake.input.llocal_bdg
d_estimate = snakemake.input.d_estimate
control_bam = snakemake.input.ctrl_bam
output_bdg = snakemake.output.bdg

if d_bdg is None:
    raise ValueError("macs2/combine_backgrounds requires input.d_bdg")
if slocal_bdg is None:
    raise ValueError("macs2/combine_backgrounds requires input.slocal_bdg")
if llocal_bdg is None:
    raise ValueError("macs2/combine_backgrounds requires input.llocal_bdg")
if d_estimate is None:
    raise ValueError("macs2/combine_backgrounds requires input.d_estimate")
else:
    with open(d_estimate) as f:
        d_estimate = f.readline().strip()
if control_bam is None:
    raise ValueError("macs2/combine_backgrounds requires input.ctrl_bam")
else:
    control_bam = subprocess.run(['samtools', 'view', '-c', control_bam],
                                 stdout=subprocess.PIPE).stdout.decode('utf-8')
if output_bdg is None:
    raise ValueError("macs2/combine_backgrounds requires output.bdg")

effective_genome_count = snakemake.params.block.get('effective_genome_count',
                        snakemake.params.block.get('reference_effective_genome_count', ''))

genome_background = float(control_bam) * float(d_estimate) / float(effective_genome_count)

tmp1 = tempfile.TemporaryFile()
tmp2 = tempfile.TemporaryFile()
cmds = (
    'macs2 bdgcmp -m max '
    '-t {slocal_bdg} '
    '-c {llocal_bdg} '
    '-o {tmp1} && '
    'macs2 bdgcmp -m max '
    '-t {tmp1} '
    '-c {d_bdg} '
    '-o {tmp2} && '
    'macs2 bdgopt -m max '
    '-i {tmp2} '
    '-m max '
    '-p {genome_background} '
    '-o {output_bdg} '
)
shell(cmds + '{log}')
