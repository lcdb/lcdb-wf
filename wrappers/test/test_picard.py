import pytest
import os
import gzip
from utils import run, dpath, rm, symlink_in_tempdir


@pytest.fixture(scope='session')
def sample1_se_bam_markdups(sample1_se_bam, tmpdir_factory):
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
            sample1_se_bam: 'sample1.bam',
        }
    )
    tmpdir = str(tmpdir_factory.mktemp('markduplicates_fixture'))
    run(dpath('../wrappers/picard/markduplicates'), snakefile, None, input_data_func, tmpdir, use_conda=True)
    return {
            'bam': os.path.join(tmpdir, 'sample1.dupsmarked.bam'),
            'metrics': os.path.join(tmpdir, 'sample1.dupmetrics.txt')
            }


def test_markduplicates_se(sample1_se_bam_markdups, tmpdir):
    assert open(sample1_se_bam_markdups['metrics']).readline().startswith('##')


def test_picard_collectrnaseqmetrics_se(sample1_se_tiny_bam, annotation_refflat, tmpdir):
    snakefile = '''
    rule collectrnaseqmetrics:
        input:
            bam='sample1.bam',
            refflat='dm6.refflat',
        output:
            metrics='sample1.metrics'
        log: 'log'
        params:
            extra="STRAND=NONE",
            java_args='-Xmx512m'
        wrapper: 'file:wrapper'
    '''
    input_data_func=symlink_in_tempdir(
        {
            sample1_se_tiny_bam: 'sample1.bam',
            annotation_refflat: 'dm6.refflat',
        }
    )

    def check():
        assert '## METRICS CLASS' in open('sample1.metrics').read()

    run(dpath('../wrappers/picard/collectrnaseqmetrics'), snakefile, check, input_data_func, tmpdir, use_conda=True)


def test_picard_collectrnaseqmetrics_se_plot(sample1_se_tiny_bam, annotation_refflat, tmpdir):
    snakefile = '''
    rule collectrnaseqmetrics:
        input:
            bam='sample1.bam',
            refflat='dm6.refflat',
        output:
            metrics='sample1.metrics',
            plot='sample1.pdf'
        log: 'log'
        params: extra="STRAND=NONE CHART=sample1.pdf"
        wrapper: 'file:wrapper'
    '''
    input_data_func=symlink_in_tempdir(
        {
            sample1_se_tiny_bam: 'sample1.bam',
            annotation_refflat: 'dm6.refflat',
        }
    )

    def check():
        assert '## METRICS CLASS' in open('sample1.metrics').read()

    run(dpath('../wrappers/picard/collectrnaseqmetrics'), snakefile, check, input_data_func, tmpdir, use_conda=True)


@pytest.mark.xfail
def test_picard_collectrnaseqmetrics_too_small_heap(sample1_se_tiny_bam, annotation_refflat, tmpdir):
    # set the java vm heap size to 128 bytes which should fail. This tests to
    # make sure the java args are making it through to the wrapper.
    snakefile = '''
    rule collectrnaseqmetrics:
        input:
            bam='sample1.bam',
            refflat='dm6.refflat',
        output:
            metrics='sample1.metrics'
        log: 'log'
        params:
            extra="STRAND=NONE",
            java_args='-Xmx128'
        wrapper: 'file:wrapper'
    '''
    input_data_func=symlink_in_tempdir(
        {
            sample1_se_tiny_bam: 'sample1.bam',
            annotation_refflat: 'dm6.refflat',
        }
    )

    def check():
        assert '## METRICS CLASS' in open('sample1.metrics').read()

    run(dpath('../wrappers/picard/collectrnaseqmetrics'), snakefile, check, input_data_func, tmpdir, use_conda=True)
