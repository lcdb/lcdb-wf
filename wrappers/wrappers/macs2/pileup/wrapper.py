import os
import os.path
from snakemake import shell
import tempfile

log = snakemake.log_fmt_shell()
logfile = None

input_bam = snakemake.input.bam
output_bdg = snakemake.output.bdg

if input_bam is None:
    raise ValueError("macs2/pileup requires input.bam")

if output_bdg is None:
    raise ValueError("macs2/pileup requires output.bdg")

cmds = (
    'macs2 pileup '
    '-i {input_bam} '
    '-f BAM '
    '-o {output_bdg} '
)

shell(cmds + ' {log}')
