import os
from snakemake.shell import shell
import sys
sys.path.append(os.path.abspath('../..'))
from lib import aligners
import tempfile

__author__ = "Ryan Dale"
__copyright__ = "Copyright 2016, Ryan Dale"
__email__ = "dalerr@niddk.nih.gov"
__license__ = "MIT"

# Pull in parameters
extra = snakemake.params.get('extra', '')
aligner = snakemake.params.get('aligner', 'bowtie2')
subset = snakemake.params.get('subset', 100000)

if aligner == 'bowtie2':
    parse_index = aligners.prefix_from_bowtie2_index

# Make log
log = snakemake.log_fmt_shell()

# snakemake.params.fastq_screen_config can be either a dict or a string. If
# string, interpret as a filename pointing to the fastq_screen config file.
# Otherwise, create a new tempfile out of the contents of the dict:

tmp = tempfile.NamedTemporaryFile(delete=False).name
with open(tmp, 'w') as fout:
    for k, v in snakemake.input.items():
        if k != 'fastq':
            label = k
            if isinstance(v, str):
                v = [v]
            index = parse_index(v)
            fout.write(
                '\t'.join(['DATABASE', label, index, aligner.upper()]) + '\n')
    config_file = tmp

# fastq_screen hard-codes filenames according to this prefix. We will send
# hard-coded output to a temp dir, and then move them later.
tempdir = tempfile.mkdtemp()

# Note that we assume only R1 is coming in.
prefix = os.path.basename(snakemake.input.fastq[0].split('.fastq')[0])

shell(
    "fastq_screen --outdir {tempdir} "
    "--force "
    "--aligner {aligner} "
    "--conf {config_file} "
    "--subset {subset} "
    "--threads {snakemake.threads} "
    "{extra} "
    "{snakemake.input.fastq} "
    "{log}"
)

# Move output to the filenames specified by the rule
shell("cp {tempdir}/{prefix}_screen.txt {snakemake.output.txt}")

# Check for the output of the --tag option to fastq_screen
if os.path.isfile("{tempdir}/{prefix}.tagged.fastq.gz"):
    shell("cp {tempdir}/{prefix}.tagged.fastq.gz {snakemake.output.txt}.tagged.fastq.gz")

# Check for the output of the --filter XXXXXX option to fastq_screen
if os.path.isfile("{tempdir}/{prefix}.tagged_filter.fastq.gz"):
    shell("cp {tempdir}/{prefix}.tagged_filter.fastq.gz {snakemake.output.txt}.tagged_filter.fastq.gz")

# Clean up temp
shell("rm -r {tempdir}")
shell("rm {tmp}")
