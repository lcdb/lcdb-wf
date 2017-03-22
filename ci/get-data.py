import os
from snakemake.shell import shell
from snakemake.utils import makedirs

URL = 'https://github.com/lcdb/lcdb-test-data/blob/master/data/{}?raw=true'

def _download_file(fn, dest=None):
    url = URL.format(fn)
    if dest is None:
        dest = fn
    makedirs(os.path.dirname(dest))
    basename = os.path.basename(fn)
    shell('wget -q -O- {url} > {dest}')
    return dest

_download_file('samples/sample1/sample1_R1.fastq.gz')
_download_file('samples/sample2/sample2_R1.fastq.gz')
