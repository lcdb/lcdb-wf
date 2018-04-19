import os
import os.path
from snakemake import shell
import subprocess

log = snakemake.log_fmt_shell()
logfile = None
extra = snakemake.params.get('extra', '')

label = snakemake.params.block['label']
extra = snakemake.params.block.get('extra', '')

effective_genome_count = snakemake.params.block.get('effective_genome_count',
                         snakemake.params.block.get('reference_effective_genome_count', ''))

mfold_lower = snakemake.params.block.get('mfold_lower', '')
mfold_upper = snakemake.params.block.get('mfold_upper', '')

genome_count_flag = ''
if effective_genome_count != '':
    genome_count_flag = ' -g ' + effective_genome_count + ' '

d_estimate_file = snakemake.output.d_estimate

# determine if this parameter set has already been estimated for this file
found_estimate = False
if os.path.isfile(d_estimate_file):
    with f as open(d_estimate_file):
        for line in f:
            parsed = line.split()
            if len(parsed) == 3 and
               parsed[0] == mfold_lower and
               parsed[1] == mfold_upper:
               found_estimate = True

if not found_estimate:
    cmds = (
        'macs2 predictd '
        '-i {snakemake.input.bam} '
        '{genome_count_flag} '
        '-m {mfold_lower} {mfold_upper} '
        '-f BAM '
    )
    shell(cmds + ' {log}')
    shell("""awk '/predicted fragment length is/ {print {0}" "{1}" "$(NF-1)}' """
          """< {2} >> {3}""".format(mfold_lower,
                                    mfold_upper,
                                    log.split(' ')[-1],
                                    d_estimate_file))
    shell("""awk '/tag size is determined as/ {print $(NF-1)}' """
          """< {0} > {1}""".format(log.split(' ')[-1],
                                   snakemake.output.frag_length))
