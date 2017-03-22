import subprocess as sp
import pytest
from snakemake import shell

def test_samtools_sort_and_index(sample1_se_sort_bam, sample1_se_sort_bam_bai):
    with pytest.raises(sp.CalledProcessError):
        shell('samtools view {sample1_se_sort_bam} 2L:1-100')
    shell('samtools view {sample1_se_sort_bam_bai[bam]} 2L:1-100')

