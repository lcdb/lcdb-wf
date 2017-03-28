import os
import zipfile
from utils import run, dpath, rm, symlink_in_tempdir
from test_bowtie2 import bowtie2_indexes

def test_fastq_screen(sample1_se_tiny_fq, bowtie2_indexes, tmpdir):
    snakefile = '''
    rule fastq_screen:
        input:
            fastq='sample1_R1.fastq.gz',
            dm6={indexes}
        output:
            txt='sample1_R1_screen.txt'
        params:
            subset=100000,
            aligner='bowtie2'
        wrapper:
            "file:wrapper"
    '''.format(indexes=bowtie2_indexes)

    input_data_func=symlink_in_tempdir(
        {
            sample1_se_tiny_fq: 'sample1_R1.fastq.gz'
        }
    )

    def check():
        with open('sample1_R1_screen.txt') as fh:
            res = fh.readlines()
            r1 = res[0].strip().split()
            r3 = res[2].strip().split()
            assert r1[-1] == '100000'
            assert r3[0] == 'dm6'


    run(dpath('../wrappers/fastq_screen'), snakefile, check, input_data_func, tmpdir)
