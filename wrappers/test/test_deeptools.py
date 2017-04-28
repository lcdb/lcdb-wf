import os
import gzip
from utils import run, dpath, rm, symlink_in_tempdir
import pyBigWig

def test_deeptools_bamCoverage(sample1_se_bam, sample1_se_bam_bai, tmpdir):
    snakefile = '''
                rule deeptools:
                    input:
                        bam='sample1.bam',
                        bai='sample1.bam.bai'
                    output: 'sample1.bw',
                    log: 'deeptools.log'
                    wrapper: "file:wrapper"
                '''
    input_data_func=symlink_in_tempdir(
        {
            sample1_se_bam: 'sample1.bam',
            sample1_se_bam_bai['bai']: 'sample1.bam.bai',
        }
    )

    def check():
        bw = pyBigWig.open('sample1.bw')
        assert bw.header()['sumData'] == 195295397
        assert bw.stats('2L')[0] == 8.242775364434165

    run(dpath('../wrappers/deeptools/bamCoverage'), snakefile, check, input_data_func, tmpdir)
