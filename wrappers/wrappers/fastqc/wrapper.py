__author__ = "Ryan Dale"
__copyright__ = "Copyright 2016, Ryan Dale"
__email__ = "dalerr@niddk.nih.gov"
__license__ = "MIT"

import os
from snakemake.shell import shell
from snakemake.utils import makedirs

# fastqc creates a zip file and an html file but the filename is hard-coded by
# replacing fastq|fastq.gz|fq|fq.gz|bam with _fastqc.zip|_fastqc.html in the
# input file's basename.
#
# So we identify that file and move it to the expected output after fastqc is
# done.

outfile = os.path.basename(snakemake.input[0])
outdir = os.path.dirname(snakemake.output.html)
if outdir == '':
    outdir = '.'

strip = ['.fastq', '.fq', '.gz', '.bam']
for s in strip:
    outfile = outfile.replace(s, '')
out_zip = os.path.join(outdir, outfile + '_fastqc.zip')
out_html = os.path.join(outdir, outfile + '_fastqc.html')

extra = snakemake.params.get('extra', '')
log = snakemake.log_fmt_shell()

shell(
    'fastqc '
    '--threads {snakemake.threads} '
    '--noextract '
    '--quiet '
    '--outdir {outdir} '
    '{extra} '
    '{snakemake.input} '
    '{log} '
)

def same_file(x, y):
    return os.path.abspath(x) == os.path.abspath(y)

if not same_file(out_zip,snakemake.output.zip):
    shell('mv {out_zip} {snakemake.output.zip}')
if not same_file(out_html, snakemake.output.html):
    shell('mv {out_html} {snakemake.output.html}')
