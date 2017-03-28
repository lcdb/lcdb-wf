import os
from snakemake.shell import shell
from snakemake.utils import makedirs

URL = 'https://github.com/lcdb/lcdb-test-data/blob/add-chipseq/data/{}?raw=true'

def _download_file(fn, dest=None):
    url = URL.format(fn)
    if dest is None:
        dest = fn
    makedirs(os.path.dirname(dest))
    basename = os.path.basename(fn)
    shell('wget -q -O- {url} > {dest}')
    return dest

_download_file('rnaseq_samples/sample1/sample1.small_R1.fastq.gz')
_download_file('rnaseq_samples/sample2/sample2.small_R1.fastq.gz')

shell('mkdir -p data/rnaseq_samples/sample{{1,2}}')
for n in [1, 2]:
    shell(
        'mv rnaseq_samples/sample{n}/sample{n}.small_R1.fastq.gz '
        'data/rnaseq_samples/sample{n}/sample{n}_R1.fastq.gz'
    )
