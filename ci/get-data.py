#!/usr/bin/env python
import os
from snakemake.shell import shell
from snakemake.utils import makedirs

shell.executable('/bin/bash')

URL = 'https://github.com/lcdb/lcdb-test-data/blob/add-chipseq/data/{}?raw=true'

def _download_file(fn, dest=None):
    url = URL.format(fn)
    if dest is None:
        dest = fn
    makedirs(os.path.dirname(dest))
    basename = os.path.basename(fn)
    shell('wget -q -O- {url} > {dest}')
    return dest

_download_file('rnaseq_samples/sample1/sample1.small_R1.fastq.gz', 'data/example_data/rnaseq_sample1.fq.gz')
_download_file('rnaseq_samples/sample2/sample2.small_R1.fastq.gz', 'data/example_data/rnaseq_sample2.fq.gz')
_download_file('rnaseq_samples/sample3/sample3.small_R1.fastq.gz', 'data/example_data/rnaseq_sample3.fq.gz')
_download_file('rnaseq_samples/sample4/sample4.small_R1.fastq.gz', 'data/example_data/rnaseq_sample4.fq.gz')
_download_file('chipseq_samples/input_1/input_1.tiny_R1.fastq.gz', 'data/example_data/chipseq_input1.fq.gz')
_download_file('chipseq_samples/ip_1/ip_1.tiny_R1.fastq.gz', 'data/example_data/chipseq_ip1.fq.gz')
_download_file('chipseq_samples/input_2/input_2.tiny_R1.fastq.gz', 'data/example_data/chipseq_input2.fq.gz')
_download_file('chipseq_samples/ip_2/ip_2.tiny_R1.fastq.gz', 'data/example_data/chipseq_ip2.fq.gz')
_download_file('chipseq_samples/ip_3/ip_3.tiny_R1.fastq.gz', 'data/example_data/chipseq_ip3.fq.gz')
_download_file('chipseq_samples/ip_4/ip_4.tiny_R1.fastq.gz', 'data/example_data/chipseq_ip4.fq.gz')
_download_file('chipseq_samples/input_3/input_3.tiny_R1.fastq.gz', 'data/example_data/chipseq_input3.fq.gz')
