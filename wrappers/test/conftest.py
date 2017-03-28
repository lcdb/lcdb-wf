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

from raw_data_fixtures import *


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
