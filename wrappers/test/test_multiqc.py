import pytest
import os
import gzip
from utils import run, dpath, rm, symlink_in_tempdir


def test_multiqc(fastqc, tmpdir):
    snakefile = '''
    rule multiqc:
        input: 'results/sample1_R1_fastqc.zip'
        output: 'multiqc.html'
        log: 'log'
        params:
            analysis_directory='results'
        wrapper: 'file:wrapper'
    '''
    input_data_func=symlink_in_tempdir(
        {
            fastqc: 'results/sample1_R1_fastqc.zip',
        }
    )

    def check():
        assert '<!DOCTYPE html>' in open('multiqc.html').readline()

    run(dpath('../wrappers/multiqc'), snakefile, check, input_data_func, tmpdir)
