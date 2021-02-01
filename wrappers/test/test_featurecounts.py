import os
import gzip
from utils import run, dpath, rm, symlink_in_tempdir

def test_featurecounts_se(sample1_se_tiny_bam, annotation, tmpdir):
    snakefile = '''
                rule featurecounts:
                    input:
                        annotation='dm6.gtf',
                        bam='sample1.bam'
                    output:
                        counts='sample1.counts',
                    log: 'featurecounts.log'
                    wrapper: "file:wrapper"
                '''
    input_data_func=symlink_in_tempdir(
        {
            sample1_se_tiny_bam: 'sample1.bam',
            annotation: 'dm6.gtf',
        }
    )

    def check():
        assert '//===================' in open('featurecounts.log').read()
        assert '# Program:featureCounts' in open('sample1.counts').readline()
        assert open('sample1.counts.summary').readline().startswith('Status')
        assert sum(1 for _ in open('sample1.counts')) == 169

    run(dpath('../wrappers/featurecounts'), snakefile, check, input_data_func, tmpdir)

def test_featurecounts_pe(sample1_pe_tiny_bam, annotation, tmpdir):
    snakefile = '''
                rule featurecounts:
                    input:
                        annotation='dm6.gtf',
                        bam='sample1.bam'
                    output:
                        counts='sample1.counts',
                    log: 'featurecounts.log'
                    params: extra='-p -P -s 1 -B --splitOnly'
                    wrapper: "file:wrapper"
                '''
    input_data_func=symlink_in_tempdir(
        {
            sample1_pe_tiny_bam: 'sample1.bam',
            annotation: 'dm6.gtf',
        }
    )

    def check():
        assert '//===================' in open('featurecounts.log').read()
        assert '# Program:featureCounts' in open('sample1.counts').readline()
        assert open('sample1.counts.summary').readline().startswith('Status')
        assert sum(1 for _ in open('sample1.counts')) == 169

        # TODO: maybe assert that below a certain level are counted when all
        # those extra arguments are used?

    run(dpath('../wrappers/featurecounts'), snakefile, check, input_data_func, tmpdir)
