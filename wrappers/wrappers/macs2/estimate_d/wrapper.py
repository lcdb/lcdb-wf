import os
import os.path
from snakemake import shell
import tempfile

log = snakemake.log_fmt_shell()
logfile = None
extra = snakemake.params.get('extra', '')

label = snakemake.params.block['label']
extra = snakemake.params.block.get('extra', '')

effective_genome_count = snakemake.params.block.get('effective_genome_count',
                         snakemake.params.block.get('reference_effective_genome_count', ''))

mfold_lower = snakemake.params.block.get('mfold_lower', '')[0]
mfold_upper = snakemake.params.block.get('mfold_upper', '')[0]

genome_count_flag = ''
if effective_genome_count != '':
    genome_count_flag = ' -g ' + effective_genome_count + ' '

d_estimate_file = snakemake.output.d_estimate

# determine if this parameter set has already been estimated for this file
found_estimate = False
if os.path.isfile(d_estimate_file):
    with open(d_estimate_file) as f:
        for line in f:
            parsed = line.split()
            if len(parsed) == 3 and \
               parsed[0] == mfold_lower and \
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
    inlog = log.split(' ')[-2]
    tmp = tempfile.NamedTemporaryFile().name
    cmds = ('awk \'/predicted fragment length is/ {{print "{mfold_lower} {mfold_upper} "$(NF-1)}}\' '
            '< {inlog} > {tmp}')
    shell(cmds)
    with open(tmp) as f:
        lines = f.readlines()
        if len(lines) < 1:
            shell("echo \"{0} {1} {2}\" >> {3}".format(mfold_lower,
                                                       mfold_upper,
                                                       300,
                                                       d_estimate_file))
        else:
            shell("echo \"{0} {1} {2}\" >> {3}".format(mfold_lower,
                                                       mfold_upper,
                                                       lines[0].strip(),
                                                       d_estimate_file))
    cmds = ("""awk '/tag size is determined as/ {{print $(NF-1)}}' """
            """< {inlog} > {snakemake.output.read_length}""")
    shell(cmds)
    shell("rm {0}".format(tmp))
