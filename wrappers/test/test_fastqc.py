import os
import zipfile
from utils import run, dpath, rm, symlink_in_tempdir

import pytest
from utils import tmpdir_for_func, _download_file

@pytest.fixture(scope='session')
def fastqc(sample1_se_tiny_fq, tmpdir_factory):
    snakefile = '''
    rule fastqc:
        input:
            fastq='sample1_R1.fastq.gz'
        output:
            html='sample1_R1_fastqc.html',
            zip='sample1_R1_fastqc.zip'
        wrapper: "file:wrapper"'''
    input_data_func = symlink_in_tempdir(
        {
            sample1_se_tiny_fq: 'sample1_R1.fastq.gz'
        }
    )
    tmpdir = str(tmpdir_factory.mktemp('fastqc_fixture'))
    run(dpath('../wrappers/fastqc'), snakefile, None, input_data_func, tmpdir)
    return os.path.join(tmpdir, 'sample1_R1_fastqc.zip')


def test_fastqc(sample1_se_tiny_fq, tmpdir):
    snakefile = '''
    rule fastqc:
        input:
            fastq='sample1_R1.fastq.gz'
        output:
            html='results/sample1_R1.html',
            zip='sample1_R1.zip'
        wrapper: "file:wrapper"'''
    input_data_func=symlink_in_tempdir(
        {
            sample1_se_tiny_fq: 'sample1_R1.fastq.gz'
        }
    )

    def check():
        assert '<html>' in open('results/sample1_R1.html').readline()
        contents = [
            'sample1_R1_fastqc/',
            'sample1_R1_fastqc/Icons/',
            'sample1_R1_fastqc/Images/',
            'sample1_R1_fastqc/Icons/fastqc_icon.png',
            'sample1_R1_fastqc/Icons/warning.png',
            'sample1_R1_fastqc/Icons/error.png',
            'sample1_R1_fastqc/Icons/tick.png',
            'sample1_R1_fastqc/summary.txt',
            'sample1_R1_fastqc/Images/per_base_quality.png',
            'sample1_R1_fastqc/Images/per_tile_quality.png',
            'sample1_R1_fastqc/Images/per_sequence_quality.png',
            'sample1_R1_fastqc/Images/per_base_sequence_content.png',
            'sample1_R1_fastqc/Images/per_sequence_gc_content.png',
            'sample1_R1_fastqc/Images/per_base_n_content.png',
            'sample1_R1_fastqc/Images/sequence_length_distribution.png',
            'sample1_R1_fastqc/Images/duplication_levels.png',
            'sample1_R1_fastqc/Images/adapter_content.png',
            'sample1_R1_fastqc/fastqc_report.html',
            'sample1_R1_fastqc/fastqc_data.txt',
            'sample1_R1_fastqc/fastqc.fo'
        ]
        for i in zipfile.ZipFile('sample1_R1.zip').namelist():
            assert i in contents

    run(dpath('../wrappers/fastqc'), snakefile, check, input_data_func, tmpdir)
