import os
from snakemake import shell
import tempfile

log = snakemake.log_fmt_shell()
logfile = None
extra = snakemake.params.get('extra', '')

label = snakemake.params.block['label']
extra = snakemake.params.block.get('extra', '')


input_bam = snakemake.input.bam
d_estimate = snakemake.input.d_estimate
override_size = snakemake.input.override_size
output_bdg = snakemake.output.bdg

if input_bam is None:
    raise ValueError("macs2/pileup requires input.bam")

if output_bdg is None:
    raise ValueError("macs2/pileup requires output.bdg")

extsize_flag = ""
if d_estimate is not None:
    with open(d_estimate) as f:
        d_estimate = f.readline().strip()
    extsize_flag = " -B --extsize " + d_estimate + " "

cmds = (
    'macs2 pileup '
    '-i {input_bam} '
    '-f BAM '
    '{extsize_flag} '
    '-o {output_bdg} '
)

shell(cmds + ' {log}')

if override_size is not None:
    prop_adj = float(d_estimate) / float(override_size)
    tmp = tempfile.TemporaryFile()
    cmds = (
        'macs2 bdgopt '
        '-i {output_bdg} '
        '-m multiply '
        '-p {prop_adj} '
        '-o {tmp} && mv {tmp} {output_bdg} '
        ' {log}.postprocess.log'
    )
    shell(cmds)
