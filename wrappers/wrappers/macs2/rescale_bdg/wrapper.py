import os
from snakemake import shell
import tempfile

log = snakemake.log_fmt_shell()
logfile = None
extra = snakemake.params.get('extra', '')

label = snakemake.params.block['label']
extra = snakemake.params.block.get('extra', '')

ip_bam = snakemake.input.ip_bam
control_bam = snakemake.input.ctrl_bam
background_bdg = snakemake.input.background_bdg
ip_bdg = snakemake.input.ip_bdg

output_ip_bdg = snakemake.output.ip_bdg
output_lambda_bdg = snakemake.output.lambda_bdg

if ip_bam is None:
    raise ValueError("macs2/rescale_bdg requires input.ip_bam")
if control_bam is None:
    raise ValueError("macs2/rescale_bdg requires input.control_bam")
if background_bdg is None:
    raise ValueError("macs2/rescale_bdg requires input.background_bdg")
if ip_bdg is None:
    raise ValueError("macs2/rescale_bdg requires input.ip_bdg")
if output_ip_bdg is None:
    raise ValueError("macs2/rescale_bdg requires output.ip_bdg")
if output_lambda_bdg is None:
    raise ValueError("macs2/rescale_bdg requires output.lambda_bdg")

control_bam = float(subprocess.run(['samtools', 'view', '-c', control_bam],
                                   stdout=subprocess.PIPE).stdout.decode('utf-8'))
ip_bam = float(subprocess.run(['samtools', 'view', '-c', ip_bam],
                              stdout=subprocess.PIPE).stdout.decode('utf-8'))

if control_bam > ip_bam:
    mult_factor = ip_bam / control_bam
    shell('ln -s {0} {1}'.format(ip_bdg,
                                 output_ip_bdg))
    cmds = (
        'macs2 bdgopt '
        '-i {background_bdg} '
        '-m multiply '
        '-p {mult_factor} '
        '-o {output_lambda_bdg} ')
else:
    mult_factor = control_bam / ip_bam
    shell('ln -s {0} {1}'.format(background_bdg,
                                 output_lamdba_bdg))

    cmds = (
        'macs2 bdgopt '
        '-i {ip_bdg} '
        '-m multiply '
        '-p {mult_factor} '
        '-o {output_ip_bdg} ')
shell(cmds + '{log}')
