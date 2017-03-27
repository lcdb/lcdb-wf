import os
import pytest
from snakemake.shell import shell
from lcdblib.snakemake import aligners
from utils import run, dpath, symlink_in_tempdir, tmpdir_for_func


@pytest.fixture(scope='session')
def bowtie2_indexes(dm6_fa, tmpdir_factory):
    d = tmpdir_for_func(tmpdir_factory)
    snakefile = '''
    rule bowtie2:
        input: fasta='dm6.fa'
        output: index=['dm6.1.bt2', 'dm6.2.bt2']
        log: 'bowtie2.log'
        wrapper: 'file:wrapper'
    '''
    input_data_func = symlink_in_tempdir(
        {
            dm6_fa: 'dm6.fa'
        }
    )

    def check():
        assert 'Total time for backward call to driver' in open('bowtie2.log').readlines()[-1]
        assert list(shell('bowtie2-inspect dm6 -n', iterable=True)) == ['2L', '2R']

    run(
        dpath('../wrappers/bowtie2/build'),
        snakefile, check, input_data_func, d)
    return aligners.bowtie2_index_from_prefix(os.path.join(d, 'dm6'))


def _dict_of_bowtie2_indexes(bowtie2_indexes, prefix):
    d = {}
    indexes = aligners.bowtie2_index_from_prefix(prefix)
    bowtie2_indexes = sorted(bowtie2_indexes)
    indexes = sorted(indexes)
    for k, v in zip(bowtie2_indexes, indexes):
        d[k] = v
    return d


def test_bowtie2_align_se(bowtie2_indexes, sample1_se_tiny_fq, tmpdir):
    d = _dict_of_bowtie2_indexes(bowtie2_indexes, 'dm6')
    indexes = list(d.values())
    snakefile = '''
        rule bowtie2_align:
            input:
                fastq='sample1_R1.fastq.gz',
                index={indexes}
            output:
                bam='sample1.bam'
            log: "bowtie2.log"
            wrapper: "file:wrapper"
    '''.format(indexes=indexes)
    d[sample1_se_tiny_fq] = 'sample1_R1.fastq.gz'
    input_data_func = symlink_in_tempdir(d)

    def check():
        assert "overall alignment rate" in open('bowtie2.log').read()

        # should have at least some mapped and unmapped
        assert int(list(shell('samtools view -c -f 0x04 sample1.bam', iterable=True))[0]) > 0
        assert int(list(shell('samtools view -c -F 0x04 sample1.bam', iterable=True))[0]) > 0

    run(dpath('../wrappers/bowtie2/align'), snakefile, check, input_data_func, tmpdir)


def test_bowtie2_align_se_rm_unmapped(bowtie2_indexes, sample1_se_tiny_fq, tmpdir):
    d = _dict_of_bowtie2_indexes(bowtie2_indexes, 'dm6')
    indexes = list(d.values())
    snakefile = '''
        rule bowtie2_align:
            input:
                fastq='sample1_R1.fastq.gz',
                index={indexes}
            output:
                bam='sample1.bam'
            params:
                samtools_view_extra='-F 0x04'
            log: "bowtie2.log"
            wrapper: "file:wrapper"
    '''.format(indexes=indexes)
    d[sample1_se_tiny_fq] = 'sample1_R1.fastq.gz'
    input_data_func = symlink_in_tempdir(d)

    def check():
        assert "overall alignment rate" in open('bowtie2.log').read()

        # should have at least some mapped and unmapped
        assert int(list(shell('samtools view -c -f 0x04 sample1.bam', iterable=True))[0]) == 0
        assert int(list(shell('samtools view -c -F 0x04 sample1.bam', iterable=True))[0]) > 0

    run(dpath('../wrappers/bowtie2/align'), snakefile, check, input_data_func, tmpdir)
