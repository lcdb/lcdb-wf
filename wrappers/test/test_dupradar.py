import os
import pytest
from test_picard import sample1_se_bam_markdups
from utils import symlink_in_tempdir, run, dpath


@pytest.fixture(scope='session')
def sample1_se_dupradar(sample1_se_bam_markdups, annotation, tmpdir_factory):
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
            dataframe='sample1.dupradar.tsv',
            model='sample1.model.txt',
            curve='sample1.curve.txt'
        wrapper:
            'file:wrapper'
    '''
    input_data_func = symlink_in_tempdir(
        {
            sample1_se_bam_markdups['bam']: 'sample1.bam',
            annotation: 'dm6.gtf',
        }
    )
    tmpdir = str(tmpdir_factory.mktemp('dupradar_fixture'))
    run(dpath('../wrappers/dupradar'), snakefile, None, input_data_func, tmpdir, use_conda=False)
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


#@pytest.mark.xfail
def test_dupradar(sample1_se_dupradar):
    assert open(sample1_se_dupradar['dataframe']).readline().startswith('"ID"\t"geneLength"')
