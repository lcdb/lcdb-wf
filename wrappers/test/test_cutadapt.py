import os
import gzip
from utils import run, dpath, rm, symlink_in_tempdir

def test_cutadapt_simple(sample1_se_tiny_fq, tmpdir):
    snakefile = '''
                rule cutadapt:
                    input:
                        fastq='sample1_R1.fastq.gz'
                    output:
                        fastq='sample1_R1.trim.fastq.gz'
                    params: extra='-a AAA'
                    wrapper: "file:wrapper"
                '''
    input_data_func=symlink_in_tempdir(
        {
            sample1_se_tiny_fq: 'sample1_R1.fastq.gz'
        }
    )

    def check():
        """
        check for line lengths and that they are at least different sized
        """
        a = sum(1 for _ in gzip.open('sample1_R1.fastq.gz'))
        b = sum(1 for _ in gzip.open('sample1_R1.trim.fastq.gz'))
        assert a == b == 4040

        assert os.path.getsize('sample1_R1.fastq.gz') != os.path.getsize('sample1_R1.trim.fastq.gz')

    run(dpath('../wrappers/cutadapt'), snakefile, check, input_data_func, tmpdir)


def test_cutadapt_simple_with_log(sample1_se_tiny_fq, tmpdir):
    snakefile = '''
                rule cutadapt:
                    input:
                        fastq='sample1_R1.fastq.gz'
                    output:
                        fastq='sample1_R1.trim.fastq.gz'
                    params: extra='-a AAA'
                    log: 'sample1.cutadapt.log'
                    wrapper: "file:wrapper"
                '''
    input_data_func=symlink_in_tempdir(
        {
            sample1_se_tiny_fq: 'sample1_R1.fastq.gz'
        }
    )

    def check():
        """
        check for line lengths and that they are at least different sized
        """
        a = sum(1 for _ in gzip.open('sample1_R1.fastq.gz'))
        b = sum(1 for _ in gzip.open('sample1_R1.trim.fastq.gz'))
        assert a == b == 4040
        assert 'This is cutadapt' in open('sample1.cutadapt.log').readline()

        assert os.path.getsize('sample1_R1.fastq.gz') != os.path.getsize('sample1_R1.trim.fastq.gz')

    run(dpath('../wrappers/cutadapt'), snakefile, check, input_data_func, tmpdir)


def test_cutadapt_se_with_list(sample1_se_tiny_fq, tmpdir):
    snakefile = '''
                rule cutadapt:
                    input: 'sample1_R1.fastq.gz'
                    output: 'sample1_R1.trim.fastq.gz'
                    params: extra='-a AAA'
                    wrapper: "file:wrapper"
                '''
    input_data_func=symlink_in_tempdir(
        {
            sample1_se_tiny_fq: 'sample1_R1.fastq.gz'
        }
    )

    def check():
        """
        check for line lengths and that they are at least different sized
        """
        a = sum(1 for _ in gzip.open('sample1_R1.fastq.gz'))
        b = sum(1 for _ in gzip.open('sample1_R1.trim.fastq.gz'))
        assert a == b == 4040

        assert os.path.getsize('sample1_R1.fastq.gz') != os.path.getsize('sample1_R1.trim.fastq.gz')

    run(dpath('../wrappers/cutadapt'), snakefile, check, input_data_func, tmpdir)

def test_cutadapt_pe(sample1_pe_tiny_fq, tmpdir):
    snakefile = '''
                rule cutadapt:
                    input:
                        R1='sample1_R1.fastq.gz',
                        R2='sample1_R2.fastq.gz',
                    output:
                        R1='sample1_R1.trim.fastq.gz',
                        R2='sample2_R1.trim.fastq.gz',
                    params: extra='-a AAA'
                    log: 'sample1.cutadapt.log'
                    wrapper: "file:wrapper"
                '''
    input_data_func=symlink_in_tempdir(
        {
            sample1_pe_tiny_fq[0]: 'sample1_R1.fastq.gz',
            sample1_pe_tiny_fq[1]: 'sample1_R2.fastq.gz',
        }
    )

    def check():
        """
        check for line lengths and that they are at least different sized
        """
        a = sum(1 for _ in gzip.open('sample1_R1.fastq.gz'))
        b = sum(1 for _ in gzip.open('sample1_R1.trim.fastq.gz'))
        assert a == b == 4040
        assert 'This is cutadapt' in open('sample1.cutadapt.log').readline()

        assert os.path.getsize('sample1_R1.fastq.gz') != os.path.getsize('sample1_R1.trim.fastq.gz')

    run(dpath('../wrappers/cutadapt'), snakefile, check, input_data_func, tmpdir)

def test_cutadapt_pe_with_list(sample1_pe_tiny_fq, tmpdir):
    snakefile = '''
                rule cutadapt:
                    input: 'sample1_R1.fastq.gz', 'sample1_R2.fastq.gz',
                    output: 'sample1_R1.trim.fastq.gz', 'sample2_R1.trim.fastq.gz',
                    params: extra='-a AAA'
                    log: 'sample1.cutadapt.log'
                    wrapper: "file:wrapper"
                '''
    input_data_func=symlink_in_tempdir(
        {
            sample1_pe_tiny_fq[0]: 'sample1_R1.fastq.gz',
            sample1_pe_tiny_fq[1]: 'sample1_R2.fastq.gz',
        }
    )

    def check():
        """
        check for line lengths and that they are at least different sized
        """
        a = sum(1 for _ in gzip.open('sample1_R1.fastq.gz'))
        b = sum(1 for _ in gzip.open('sample1_R1.trim.fastq.gz'))
        assert a == b == 4040
        assert 'This is cutadapt' in open('sample1.cutadapt.log').readline()

        assert os.path.getsize('sample1_R1.fastq.gz') != os.path.getsize('sample1_R1.trim.fastq.gz')

    run(dpath('../wrappers/cutadapt'), snakefile, check, input_data_func, tmpdir)
