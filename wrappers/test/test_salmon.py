import os
import pytest
from snakemake.shell import shell
from utils import run, dpath, rm, symlink_in_tempdir, tmpdir_for_func


@pytest.fixture(scope='session')
def salmon_index(tmpdir_factory, transcriptome):
    d = tmpdir_for_func(tmpdir_factory)
    snakefile = '''
    rule salmon:
        input: fasta='transcriptome.fa'
        output: hash='salmon_index/hash.bin'
        log: 'log'
        wrapper: 'file:wrapper'
    '''
    input_data_func = symlink_in_tempdir(
        {
            transcriptome: 'transcriptome.fa',
        }
    )

    def check():
        log = open('log').read()
        assert '[info] done building index' in log

    run(
        dpath('../wrappers/salmon/index'),
        snakefile, check, input_data_func, d)
    return os.path.join(d, 'salmon_index')


def test_salmon_quant(tmpdir, sample1_se_tiny_fq, salmon_index):
    snakefile = '''
    rule salmon_quant:
        input:
             unmatedReads='sample1.fq.gz',
             index=['idx/hash.bin', 'idx/sa.bin']
        output: 'sample1/salmon/quant.sf'
        params: extra='--libType A'
        log: 'salmon.quant.log'
        wrapper: 'file:wrapper'
    '''
    input_data_func = symlink_in_tempdir(
        {
            sample1_se_tiny_fq: 'sample1.fq.gz',
            salmon_index: 'idx',
        }
    )

    def check():
        assert open('sample1/salmon/quant.sf').readline() == (
                'Name\tLength\tEffectiveLength\tTPM\tNumReads\n')

    run(
        dpath('../wrappers/salmon/quant'),
        snakefile, check, input_data_func, tmpdir)

def test_salmon_quant_single_index(tmpdir, sample1_se_tiny_fq, salmon_index):
    snakefile = '''
    rule salmon_quant:
        input:
             unmatedReads='sample1.fq.gz',
             index='idx/hash.bin'
        output: 'sample1/salmon/quant.sf'
        params: extra='--libType A'
        log: 'salmon.quant.log'
        wrapper: 'file:wrapper'
    '''
    input_data_func = symlink_in_tempdir(
        {
            sample1_se_tiny_fq: 'sample1.fq.gz',
            salmon_index: 'idx',
        }
    )

    def check():
        assert open('sample1/salmon/quant.sf').readline() == (
                'Name\tLength\tEffectiveLength\tTPM\tNumReads\n')

    run(
        dpath('../wrappers/salmon/quant'),
        snakefile, check, input_data_func, tmpdir)
