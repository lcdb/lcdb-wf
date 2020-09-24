#!/usr/bin/env python
import os
from snakemake.shell import shell
from snakemake.utils import makedirs

shell.executable('/bin/bash')
BRANCH = 'master'
URL = 'https://github.com/lcdb/lcdb-test-data-human/blob/{0}/{{}}?raw=true'.format(BRANCH)


def _download_file(fn, dest=None):
    url = URL.format(fn)
    if dest is None:
        dest = fn
    makedirs(os.path.dirname(dest))
    basename = os.path.basename(fn)
    shell('wget -q -O- {url} > {dest}')
    return dest

_download_file('rnaseq_samples/A549_DMSO_1h_1/A549_DMSO_1h_1.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/A549_DMSO_1h_1.small_R1.fastq.gz')
_download_file('rnaseq_samples/A549_DMSO_1h_2/A549_DMSO_1h_2.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/A549_DMSO_1h_2.small_R1.fastq.gz')
_download_file('rnaseq_samples/A549_DMSO_6h_1/A549_DMSO_6h_1.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/A549_DMSO_6h_1.small_R1.fastq.gz')
_download_file('rnaseq_samples/A549_DMSO_6h_2/A549_DMSO_6h_2.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/A549_DMSO_6h_2.small_R1.fastq.gz')
_download_file('rnaseq_samples/A549_JQ1_1h_1/A549_JQ1_1h_1.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/A549_JQ1_1h_1.small_R1.fastq.gz')
_download_file('rnaseq_samples/A549_JQ1_1h_2/A549_JQ1_1h_2.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/A549_JQ1_1h_2.small_R1.fastq.gz')
_download_file('rnaseq_samples/A549_JQ1_6h_1/A549_JQ1_6h_1.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/A549_JQ1_6h_1.small_R1.fastq.gz')
_download_file('rnaseq_samples/A549_JQ1_6h_2/A549_JQ1_6h_2.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/A549_JQ1_6h_2.small_R1.fastq.gz')
_download_file('rnaseq_samples/A549_dBet6_1h_1/A549_dBet6_1h_1.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/A549_dBet6_1h_1.small_R1.fastq.gz')
_download_file('rnaseq_samples/A549_dBet6_1h_2/A549_dBet6_1h_2.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/A549_dBet6_1h_2.small_R1.fastq.gz')
_download_file('rnaseq_samples/A549_dBet6_6h_1/A549_dBet6_6h_1.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/A549_dBet6_6h_1.small_R1.fastq.gz')
_download_file('rnaseq_samples/A549_dBet6_6h_2/A549_dBet6_6h_2.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/A549_dBet6_6h_2.small_R1.fastq.gz')
_download_file('rnaseq_samples/HAP1_DMSO_1h_1/HAP1_DMSO_1h_1.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/HAP1_DMSO_1h_1.small_R1.fastq.gz')
_download_file('rnaseq_samples/HAP1_DMSO_1h_2/HAP1_DMSO_1h_2.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/HAP1_DMSO_1h_2.small_R1.fastq.gz')
_download_file('rnaseq_samples/HAP1_DMSO_6h_1/HAP1_DMSO_6h_1.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/HAP1_DMSO_6h_1.small_R1.fastq.gz')
_download_file('rnaseq_samples/HAP1_DMSO_6h_2/HAP1_DMSO_6h_2.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/HAP1_DMSO_6h_2.small_R1.fastq.gz')
_download_file('rnaseq_samples/HAP1_JQ1_1h_1/HAP1_JQ1_1h_1.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/HAP1_JQ1_1h_1.small_R1.fastq.gz')
_download_file('rnaseq_samples/HAP1_JQ1_1h_2/HAP1_JQ1_1h_2.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/HAP1_JQ1_1h_2.small_R1.fastq.gz')
_download_file('rnaseq_samples/HAP1_JQ1_6h_1/HAP1_JQ1_6h_1.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/HAP1_JQ1_6h_1.small_R1.fastq.gz')
_download_file('rnaseq_samples/HAP1_JQ1_6h_2/HAP1_JQ1_6h_2.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/HAP1_JQ1_6h_2.small_R1.fastq.gz')
_download_file('rnaseq_samples/HAP1_dBet6_1h_1/HAP1_dBet6_1h_1.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/HAP1_dBet6_1h_1.small_R1.fastq.gz')
_download_file('rnaseq_samples/HAP1_dBet6_1h_2/HAP1_dBet6_1h_2.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/HAP1_dBet6_1h_2.small_R1.fastq.gz')
_download_file('rnaseq_samples/HAP1_dBet6_6h_1/HAP1_dBet6_6h_1.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/HAP1_dBet6_6h_1.small_R1.fastq.gz')
_download_file('rnaseq_samples/HAP1_dBet6_6h_2/HAP1_dBet6_6h_2.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/HAP1_dBet6_6h_2.small_R1.fastq.gz')
_download_file('rnaseq_samples/K562_DMSO_1h_1/K562_DMSO_1h_1.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/K562_DMSO_1h_1.small_R1.fastq.gz')
_download_file('rnaseq_samples/K562_DMSO_1h_2/K562_DMSO_1h_2.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/K562_DMSO_1h_2.small_R1.fastq.gz')
_download_file('rnaseq_samples/K562_DMSO_6h_1/K562_DMSO_6h_1.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/K562_DMSO_6h_1.small_R1.fastq.gz')
_download_file('rnaseq_samples/K562_DMSO_6h_2/K562_DMSO_6h_2.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/K562_DMSO_6h_2.small_R1.fastq.gz')
_download_file('rnaseq_samples/K562_JQ1_1h_1/K562_JQ1_1h_1.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/K562_JQ1_1h_1.small_R1.fastq.gz')
_download_file('rnaseq_samples/K562_JQ1_1h_2/K562_JQ1_1h_2.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/K562_JQ1_1h_2.small_R1.fastq.gz')
_download_file('rnaseq_samples/K562_JQ1_6h_1/K562_JQ1_6h_1.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/K562_JQ1_6h_1.small_R1.fastq.gz')
_download_file('rnaseq_samples/K562_JQ1_6h_2/K562_JQ1_6h_2.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/K562_JQ1_6h_2.small_R1.fastq.gz')
_download_file('rnaseq_samples/K562_dBet6_1h_2/K562_dBet6_1h_2.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/K562_dBet6_1h_2.small_R1.fastq.gz')
_download_file('rnaseq_samples/K562_dBet6_6h_1/K562_dBet6_6h_1.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/K562_dBet6_6h_1.small_R1.fastq.gz')
_download_file('rnaseq_samples/K562_dBet6_6h_2/K562_dBet6_6h_2.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/K562_dBet6_6h_2.small_R1.fastq.gz')
_download_file('rnaseq_samples/outlier_HAP1_dBet6_1h_1/outlier_HAP1_dBet6_1h_1.small_R1.fastq.gz', 'workflows/rnaseq/data/example_data/outlier_HAP1_dBet6_1h_1.small_R1.fastq.gz')

_download_file('chipseq_samples/BRD4_DMSO_1/BRD4_DMSO_1.small_R1.fastq.gz', 'workflows/chipseq/data/example_data/BRD4_DMSO_1.small_R1.fastq.gz')
_download_file('chipseq_samples/BRD4_DMSO_2/BRD4_DMSO_2.small_R1.fastq.gz', 'workflows/chipseq/data/example_data/BRD4_DMSO_2.small_R1.fastq.gz')
_download_file('chipseq_samples/BRD4_dBET6_1/BRD4_dBET6_1.small_R1.fastq.gz', 'workflows/chipseq/data/example_data/BRD4_dBET6_1.small_R1.fastq.gz')
_download_file('chipseq_samples/BRD4_dBET6_2/BRD4_dBET6_2.small_R1.fastq.gz', 'workflows/chipseq/data/example_data/BRD4_dBET6_2.small_R1.fastq.gz')
_download_file('chipseq_samples/MTHFD1_DMSO_1/MTHFD1_DMSO_1.small_R1.fastq.gz', 'workflows/chipseq/data/example_data/MTHFD1_DMSO_1.small_R1.fastq.gz')
_download_file('chipseq_samples/MTHFD1_DMSO_2/MTHFD1_DMSO_2.small_R1.fastq.gz', 'workflows/chipseq/data/example_data/MTHFD1_DMSO_2.small_R1.fastq.gz')
_download_file('chipseq_samples/MTHFD1_dBET6_1/MTHFD1_dBET6_1.small_R1.fastq.gz', 'workflows/chipseq/data/example_data/MTHFD1_dBET6_1.small_R1.fastq.gz')
_download_file('chipseq_samples/MTHFD1_dBET6_2/MTHFD1_dBET6_2.small_R1.fastq.gz', 'workflows/chipseq/data/example_data/MTHFD1_dBET6_2.small_R1.fastq.gz')
_download_file('chipseq_samples/input_DMSO_1/input_DMSO_1.small_R1.fastq.gz', 'workflows/chipseq/data/example_data/input_DMSO_1.small_R1.fastq.gz')
_download_file('chipseq_samples/input_DMSO_2/input_DMSO_2.small_R1.fastq.gz', 'workflows/chipseq/data/example_data/input_DMSO_2.small_R1.fastq.gz')
_download_file('chipseq_samples/input_dBET6_1/input_dBET6_1.small_R1.fastq.gz', 'workflows/chipseq/data/example_data/input_dBET6_1.small_R1.fastq.gz')
_download_file('chipseq_samples/input_dBET6_2/input_dBET6_2.small_R1.fastq.gz', 'workflows/chipseq/data/example_data/input_dBET6_2.small_R1.fastq.gz')
_download_file('chipseq_samples/mockIgG_DMSO_1/mockIgG_DMSO_1.small_R1.fastq.gz', 'workflows/chipseq/data/example_data/mockIgG_DMSO_1.small_R1.fastq.gz')
_download_file('chipseq_samples/mockIgG_DMSO_2/mockIgG_DMSO_2.small_R1.fastq.gz', 'workflows/chipseq/data/example_data/mockIgG_DMSO_2.small_R1.fastq.gz')
_download_file('chipseq_samples/mockIgG_dBET6_1/mockIgG_dBET6_1.small_R1.fastq.gz', 'workflows/chipseq/data/example_data/mockIgG_dBET6_1.small_R1.fastq.gz')
_download_file('chipseq_samples/mockIgG_dBET6_2/mockIgG_dBET6_2.small_R1.fastq.gz', 'workflows/chipseq/data/example_data/mockIgG_dBET6_2.small_R1.fastq.gz')

