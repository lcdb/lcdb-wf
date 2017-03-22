import pytest
import pysam
from snakemake.shell import shell
from lcdblib.snakemake import aligners
from utils import run, dpath, rm, symlink_in_tempdir


def test_bowtie2_build(dm6_fa, tmpdir):
    snakefile = '''
                rule bowtie2_build:
                    input:
                        fasta='2L.fa'
                    output:
                        index=expand('data/assembly/assembly.{n}.bt2', n=range(1,5))
                    log: 'bowtie2.log'
                    wrapper: "file:wrapper"

                '''
    input_data_func=symlink_in_tempdir(
        {
            dm6_fa: '2L.fa'
        }
    )

    def check():
        assert 'Total time for backward call to driver' in open('bowtie2.log').readlines()[-1]
        assert list(shell('bowtie2-inspect data/assembly/assembly -n', iterable=True)) == ['2L']

    run(dpath('../wrappers/bowtie2/build'), snakefile, check, input_data_func, tmpdir)


def _dict_of_bowtie2_indexes(bowtie2_indexes, prefix):
    d = {}
    indexes = aligners.bowtie2_index_from_prefix(prefix)
    bowtie2_indexes = sorted(bowtie2_indexes)
    indexes = sorted(indexes)
    for k, v in zip(bowtie2_indexes, indexes):
        d[k] = v
    return d


def test_bowtie2_align_se(bowtie2_indexes, sample1_se_fq, tmpdir):
    d = _dict_of_bowtie2_indexes(bowtie2_indexes, '2L')
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
    d[sample1_se_fq] = 'sample1_R1.fastq.gz'
    input_data_func = symlink_in_tempdir(d)

    def check():
        assert "overall alignment rate" in open('bowtie2.log').read()

        # should have at least some mapped and unmapped
        assert int(list(shell('samtools view -c -f 0x04 sample1.bam', iterable=True))[0]) > 0
        assert int(list(shell('samtools view -c -F 0x04 sample1.bam', iterable=True))[0]) > 0

    run(dpath('../wrappers/bowtie2/align'), snakefile, check, input_data_func, tmpdir)


def test_bowtie2_align_se_rm_unmapped(bowtie2_indexes, sample1_se_fq, tmpdir):
    d = _dict_of_bowtie2_indexes(bowtie2_indexes, '2L')
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
    d[sample1_se_fq] = 'sample1_R1.fastq.gz'
    input_data_func = symlink_in_tempdir(d)

    def check():
        assert "overall alignment rate" in open('bowtie2.log').read()

        # should have at least some mapped and unmapped
        assert int(list(shell('samtools view -c -f 0x04 sample1.bam', iterable=True))[0]) == 0
        assert int(list(shell('samtools view -c -F 0x04 sample1.bam', iterable=True))[0]) > 0

    run(dpath('../wrappers/bowtie2/align'), snakefile, check, input_data_func, tmpdir)
