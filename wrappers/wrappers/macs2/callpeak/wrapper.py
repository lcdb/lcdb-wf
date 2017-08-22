import os
import glob
from snakemake import shell

log = snakemake.log_fmt_shell()
extra = snakemake.params.get('extra', '')

outdir, basebed = os.path.split(snakemake.output.bed)
label = snakemake.params.block['label']
extra = snakemake.params.block.get('extra', '')

cmds = (
    'macs2 '
    'callpeak '
    '-c {snakemake.input.control} '
    '-t {snakemake.input.ip} '
    '-f BAM '
    '--outdir {outdir} '
    '--name {label} '
)
# add any per-peak-calling-run extra commands
cmds += extra

shell(cmds + ' {log}')

# Depending on whether --broad was used, we may have a "*_peaks.narrowPeak" or
# a "*_peaks.broadPeak". Figure out which it is, and symlink to the output.
hits = glob.glob(os.path.join(outdir, '*_peaks.*Peak'))
assert (len(hits) == 1) and (hits[0].split('.')[-1].lower() in ['narrowpeak', 'broadpeak']), hits
basehit = os.path.basename(hits[0])
shell('cd {outdir} && ln -sf {basehit} {basebed}')
