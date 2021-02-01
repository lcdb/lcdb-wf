"""
Fixtures used for downloading data from the test data repo
"""

import os
import pytest
from utils import tmpdir_for_func, _download_file, symlink_in_tempdir, run, dpath

# ----------------------------------------------------------------------------
# FASTQ files
@pytest.fixture(scope='session')
def sample1_se_fq(tmpdir_factory):
    d = tmpdir_for_func(tmpdir_factory)
    fn = 'rnaseq_samples/sample1/sample1.small_R1.fastq.gz'
    return _download_file(fn, d)

@pytest.fixture(scope='session')
def sample1_se_tiny_fq(tmpdir_factory):
    """
    Single-end FASTQ file with 1010 reads
    """
    d = tmpdir_for_func(tmpdir_factory)
    fn = 'rnaseq_samples/sample1/sample1.tiny_R1.fastq.gz'
    return _download_file(fn, d)

@pytest.fixture(scope='session')
def sample1_pe_fq(tmpdir_factory):
    pair = []
    d = tmpdir_for_func(tmpdir_factory)
    for fn in [
        'rnaseq_samples/sample1/sample1.small_R1.fastq.gz',
        'rnaseq_samples/sample1/sample1.small_R2.fastq.gz'
    ]:
        pair.append(_download_file(fn, d))
    return pair

@pytest.fixture(scope='session')
def sample1_pe_tiny_fq(tmpdir_factory):
    pair = []
    d = tmpdir_for_func(tmpdir_factory)
    for fn in [
        'rnaseq_samples/sample1/sample1.tiny_R1.fastq.gz',
        'rnaseq_samples/sample1/sample1.tiny_R2.fastq.gz'
    ]:
        pair.append(_download_file(fn, d))
    return pair

# ----------------------------------------------------------------------------
# BAM files

@pytest.fixture(scope='session')
def sample1_se_bam(tmpdir_factory):
    d = tmpdir_for_func(tmpdir_factory)
    fn = 'rnaseq_samples/sample1/sample1.small.single.sorted.bam'
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def sample1_pe_bam(tmpdir_factory):
    d = tmpdir_for_func(tmpdir_factory)
    fn = 'rnaseq_samples/sample1/sample1.small.paired.sorted.bam'
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def sample1_se_tiny_bam(tmpdir_factory):
    d = tmpdir_for_func(tmpdir_factory)
    fn = 'rnaseq_samples/sample1/sample1.tiny.single.sorted.bam'
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def sample1_pe_tiny_bam(tmpdir_factory):
    d = tmpdir_for_func(tmpdir_factory)
    fn = 'rnaseq_samples/sample1/sample1.tiny.paired.sorted.bam'
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def sample1_se_bam_bai(sample1_se_bam, tmpdir_factory):
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
            sample1_se_bam: 'sample1.sorted.bam'

        }
    )
    tmpdir = str(tmpdir_factory.mktemp('sample1_se_bam_bai'))
    run(dpath('../wrappers/samtools/index'), snakefile, None, input_data_func, tmpdir)
    return {
            'bam': os.path.join(tmpdir, 'sample1.sorted.bam'),
            'bai': os.path.join(tmpdir, 'sample1.sorted.bam.bai'),
    }


@pytest.fixture(scope='session')
def sample1_se_tiny_bam_bai(sample1_se_tiny_bam, tmpdir_factory):
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
            sample1_se_tiny_bam: 'sample1.sorted.bam'

        }
    )
    tmpdir = str(tmpdir_factory.mktemp('sample1_se_tiny_bam_bai'))
    run(dpath('../wrappers/samtools/index'), snakefile, None, input_data_func, tmpdir)
    return {
            'bam': os.path.join(tmpdir, 'sample1.sorted.bam'),
            'bai': os.path.join(tmpdir, 'sample1.sorted.bam.bai'),
    }

# ----------------------------------------------------------------------------
# Annotations

@pytest.fixture(scope='session')
def transcriptome(tmpdir_factory):
    d = tmpdir_for_func(tmpdir_factory)
    fn = 'seq/dm6.small.transcriptome.fa'
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def dm6_fa(tmpdir_factory):
    fn = 'seq/dm6.small.fa'
    d = tmpdir_for_func(tmpdir_factory)
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def annotation(tmpdir_factory):
    fn = 'annotation/dm6.small.gtf'
    d = tmpdir_for_func(tmpdir_factory)
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def annotation_refflat(tmpdir_factory):
    fn = 'annotation/dm6.small.refflat'
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
