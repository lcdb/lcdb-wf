import os
from snakemake import shell

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

cmds = (
    'macs2 predictd '
    '-i {snakemake.input.bam} '
    '{genome_count_flag} '
    '-m {mfold_lower} {mfold_upper} '
    '-f BAM '
)

shell(cmds + ' {log}')
shell("""awk '/tag size is determined as/ {print $(NF-1)}' """
      """< {0} > {1}""".format(log.split(' ')[-1],
                               snakemake.output.frag_length))
shell("""awk '/predicted fragment length is/ {print $(NF-1)}' """
      """< {0} > {1}""".format(log.split(' ')[-1],
                               snakemake.output.d_estimate))
