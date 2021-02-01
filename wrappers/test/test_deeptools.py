import os
import gzip
from utils import run, dpath, rm, symlink_in_tempdir
import pyBigWig

def test_deeptools_bamCoverage(sample1_se_tiny_bam, sample1_se_tiny_bam_bai, tmpdir):
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
            sample1_se_tiny_bam: 'sample1.bam',
            sample1_se_tiny_bam_bai['bai']: 'sample1.bam.bai',
        }
    )

    def check():
        bw = pyBigWig.open('sample1.bw')
        header_keys = list(bw.header().keys())
        for k in ['maxVal', 'minVal', 'nBasesCovered', 'nLevels', 'sumData',
                  'sumSquared', 'version']:
            assert k in header_keys

        # bigWig version should be independent of BAM input, so we can check
        # the value
        assert bw.header()['version'] == 4

        first_chrom = list(bw.chroms().keys())[0]
        assert isinstance(bw.stats(first_chrom)[0], float)

    run(dpath('../wrappers/deeptools/bamCoverage'), snakefile, check, input_data_func, tmpdir)
