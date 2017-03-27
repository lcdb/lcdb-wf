import os
import pytest
import tempfile
import shutil
import inspect
from snakemake.shell import shell
from snakemake.utils import makedirs
from lcdblib.snakemake import aligners

from utils import (
    run,
    dpath,
    symlink_in_tempdir,
    tmpdir_for_func,
    _download_file,
)


@pytest.fixture(scope='session')
def transcriptome(tmpdir_factory):
    d = tmpdir_for_func(tmpdir_factory)
    fn = 'seq/transcriptome.fa'
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def kallisto_index(tmpdir_factory, transcriptome):
    d = tmpdir_for_func(tmpdir_factory)
    snakefile = '''
    rule kallisto:
        input: fasta='transcriptome.fa'
        output: index='transcriptome.idx'
        wrapper: 'file:wrapper'
    '''
    input_data_func = symlink_in_tempdir(
        {
            transcriptome: 'transcriptome.fa',
        }
    )

    run(
        dpath('../wrappers/kallisto/index'),
        snakefile, None, input_data_func, d)
    return os.path.join(d, 'transcriptome.idx')


@pytest.fixture(scope='session')
def sample1_se_fq(tmpdir_factory):
    d = tmpdir_for_func(tmpdir_factory)
    fn = 'samples/sample1/sample1_R1.fastq.gz'
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def sample1_pe_fq(tmpdir_factory):
    pair = []
    d = tmpdir_for_func(tmpdir_factory)
    for fn in [
        'samples/sample1/sample1_R1.fastq.gz',
        'samples/sample1/sample1_R2.fastq.gz'
    ]:
        pair.append(_download_file(fn, d))
    return pair


@pytest.fixture(scope='session')
def dm6_fa(tmpdir_factory):
    fn = 'seq/2L.fa'
    d = tmpdir_for_func(tmpdir_factory)
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def sample1_se_bam(tmpdir_factory):
    d = tmpdir_for_func(tmpdir_factory)
    fn = 'samples/sample1/sample1.single.bam'
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def sample1_se_sort_bam(sample1_se_bam, tmpdir_factory):
    snakefile = '''
    rule sort:
        input: bam='sample1.bam'
        output: bam='sample1.sorted.bam'
        wrapper: 'file:wrapper'
    '''
    input_data_func = symlink_in_tempdir(
        {
            sample1_se_bam: 'sample1.bam'

        }
    )
    tmpdir = str(tmpdir_factory.mktemp('sample1_se_sort_bam'))
    run(dpath('../wrappers/samtools/sort'), snakefile, None, input_data_func, tmpdir)
    return os.path.join(tmpdir, 'sample1.sorted.bam')


@pytest.fixture(scope='session')
def sample1_se_sort_bam_bai(sample1_se_sort_bam, tmpdir_factory):
    """
    Returns both the bam and the bam.bai
    """
    snakefile = '''
    rule index:
        input: bam='sample1.sorted.bam'
        output: bai='sample1.sorted.bam.bai'
        wrapper: 'file:wrapper'
    '''
    input_data_func = symlink_in_tempdir(
        {
            sample1_se_sort_bam: 'sample1.sorted.bam'

        }
    )
    tmpdir = str(tmpdir_factory.mktemp('sample1_se_sort_bam_bai'))
    run(dpath('../wrappers/samtools/index'), snakefile, None, input_data_func, tmpdir)
    return {
            'bam': os.path.join(tmpdir, 'sample1.sorted.bam'),
            'bai': os.path.join(tmpdir, 'sample1.sorted.bam.bai'),
    }

@pytest.fixture(scope='session')
def sample1_pe_bam(tmpdir_factory):
    d = tmpdir_for_func(tmpdir_factory)
    fn = 'samples/sample1/sample1.paired.bam'
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def sample1_pe_hisat2_bam(tmpdir_factory):
    d = tmpdir_for_func(tmpdir_factory)


@pytest.fixture(scope='session')
def hisat2_indexes(dm6_fa, tmpdir_factory):
    d = tmpdir_for_func(tmpdir_factory)
    snakefile = '''
    rule hisat2:
        input: fasta='2L.fa'
        output: index=['2L.1.ht2', '2L.2.ht2']
        wrapper: 'file:wrapper'
    '''
    input_data_func = symlink_in_tempdir(
        {
            dm6_fa: '2L.fa'
        }
    )

    run(
        dpath('../wrappers/hisat2/build'),
        snakefile, None, input_data_func, d)
    return aligners.hisat2_index_from_prefix(os.path.join(d, '2L'))


@pytest.fixture(scope='session')
def bowtie2_indexes(dm6_fa, tmpdir_factory):
    d = tmpdir_for_func(tmpdir_factory)
    snakefile = '''
    rule bowtie2:
        input: fasta='2L.fa'
        output: index=['2L.1.bt2', '2L.2.bt2']
        wrapper: 'file:wrapper'
    '''
    input_data_func = symlink_in_tempdir(
        {
            dm6_fa: '2L.fa'
        }
    )

    run(
        dpath('../wrappers/bowtie2/build'),
        snakefile, None, input_data_func, d)
    return aligners.bowtie2_index_from_prefix(os.path.join(d, '2L'))


@pytest.fixture(scope='session')
def annotation(tmpdir_factory):
    fn = 'annotation/dm6.gtf'
    d = tmpdir_for_func(tmpdir_factory)
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def annotation_refflat(tmpdir_factory):
    fn = 'annotation/dm6.refflat'
    d = tmpdir_for_func(tmpdir_factory)
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def annotation_db(annotation):
    import gffutils
    gffutils.create_db(
        data=annotation, dbfn=annotation + '.db',
        merge_strategy='merge',
        id_spec={'transcript': ['transcript_id', 'transcript_symbol'],
                 'gene': ['gene_id', 'gene_symbol']},
        gtf_transcript_key='transcript_id',
        gtf_gene_key='gene_id')

    return annotation + '.db'


@pytest.fixture(scope='session')
def annotation_bed12(annotation_db):
    import gffutils
    db = gffutils.FeatureDB(annotation_db)

    bed12 = '.'.join(annotation_db.strip().split('.')[:-2]) + '.bed12'

    with open(bed12, 'w') as handle:
        for t in db.features_of_type('transcript'):
            handle.write(db.bed12(t, name_field='transcript_id') + '\n')

    return bed12


@pytest.fixture(scope='session')
def fastqc(sample1_se_fq, tmpdir_factory):
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
            sample1_se_fq: 'sample1_R1.fastq.gz'
        }
    )
    tmpdir = str(tmpdir_factory.mktemp('fastqc_fixture'))
    run(dpath('../wrappers/fastqc'), snakefile, None, input_data_func, tmpdir)
    return os.path.join(tmpdir, 'sample1_R1_fastqc.zip')


@pytest.fixture(scope='session')
def sample1_se_bam_sorted_markdups(sample1_se_sort_bam, tmpdir_factory):
    snakefile = '''
    rule markduplicates:
        input:
            bam='sample1.bam'
        output:
            bam='sample1.dupsmarked.bam',
            metrics='sample1.dupmetrics.txt'
        log: 'log'
        wrapper: 'file:wrapper'
    '''
    input_data_func = symlink_in_tempdir(
        {
            sample1_se_sort_bam: 'sample1.bam',
        }
    )
    tmpdir = str(tmpdir_factory.mktemp('markduplicates_fixture'))
    run(dpath('../wrappers/picard/markduplicates'), snakefile, None, input_data_func, tmpdir)
    return {
            'bam': os.path.join(tmpdir, 'sample1.dupsmarked.bam'),
            'metrics': os.path.join(tmpdir, 'sample1.dupmetrics.txt')
            }


@pytest.fixture(scope='session')
def sample1_se_dupradar(sample1_se_bam_sorted_markdups, annotation, tmpdir_factory):
    snakefile = '''
    rule dupradar:
        input:
            bam='sample1.bam',
            annotation='dm6.gtf'
        output:
            density_scatter='sample1.density_scatter.png',
            expression_histogram='sample1.expression_histogram.png',
            expression_barplot='sample1.expression_barplot.png',
            expression_boxplot='sample1.expression_boxplot.png',
            multimapping_histogram='sample1.multimapping_histogram.png',
            dataframe='sample1.dupradar.tsv'
        wrapper:
            'file:wrapper'
    '''
    input_data_func = symlink_in_tempdir(
        {
            sample1_se_bam_sorted_markdups['bam']: 'sample1.bam',
            annotation: 'dm6.gtf',
        }
    )
    tmpdir = str(tmpdir_factory.mktemp('dupradar_fixture'))
    run(dpath('../wrappers/dupradar'), snakefile, None, input_data_func, tmpdir)
    mapping = dict(
        density_scatter='sample1.density_scatter.png',
        expression_histogram='sample1.expression_histogram.png',
        expression_barplot='sample1.expression_barplot.png',
        expression_boxplot='sample1.expression_boxplot.png',
        multimapping_histogram='sample1.multimapping_histogram.png',
        dataframe='sample1.dupradar.tsv',
    )
    for k, v in mapping.items():
        mapping[k] = os.path.join(tmpdir, v)
    return mapping
