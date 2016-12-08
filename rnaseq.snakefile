import os
import pandas as pd
from lcdblib.snakemake import helpers, aligners
from lcdblib.utils import utils
import common

shell.prefix('set -euo pipefail; ')

include: 'references.snakefile'

references_dir = os.environ.get('REFERENCES_DIR', config.get('references_dir', None))
if references_dir is None:
    raise ValueError('No references dir specified')
config['references_dir'] = references_dir

sampletable = pd.read_table(config['sampletable'])
samples = sampletable.ix[:, 0]

assembly = config['assembly']
refdict = common.references_dict(config)[assembly]

patterns = {
    'fastq':   'samples/{sample}/{sample}_R1.fastq.gz',
    'cutadapt': 'samples/{sample}/{sample}_R1.cutadapt.fastq.gz',
    'bam':     'samples/{sample}/{sample}.cutadapt.bam',
    'fastqc': {
        'raw': 'samples/{sample}/fastqc/{sample}_R1.fastq.gz_fastqc.zip',
        'cutadapt': 'samples/{sample}/fastqc/{sample}_R1.cutadapt.fastq.gz_fastqc.zip',
        'bam': 'samples/{sample}/fastqc/{sample}.cutadapt.bam_fastqc.zip',
    },
    'libsizes': {
        'fastq':   'samples/{sample}/{sample}_R1.fastq.gz.libsize',
        'cutadapt': 'samples/{sample}/{sample}_R1.cutadapt.fastq.gz.libsize',
        'bam':     'samples/{sample}/{sample}.cutadapt.bam.libsize',
        'rRNA':     'samples/{sample}/{sample}.cutadapt.rRNA.bam.libsize',
    },
    'rRNA': 'samples/{sample}/{sample}.cutadapt.rRNA.bam',
    'featurecounts': 'samples/{sample}/{sample}.cutadapt.bam.featurecounts.txt',
    'libsizes_table': 'libsizes_table.tsv',
    'multiqc': 'multiqc.html',
    'kallisto': {
        'h5': 'samples/{sample}/{sample}/kallisto/abundance.h5',
    },
}
fill = dict(sample=samples, count=['.count', ''])
targets = helpers.fill_patterns(patterns, fill)


def wrapper_for(path):
    return 'file://' + os.path.join('wrappers', 'wrappers', path)


rule targets:
    input:
        (
            targets['bam'] +
            utils.flatten(targets['fastqc']) +
            utils.flatten(targets['libsizes']) +
            [targets['libsizes_table']] +
            [targets['multiqc']] +
            utils.flatten(targets['featurecounts'])
        )


rule cutadapt:
    input:
        fastq=patterns['fastq']
    output:
        fastq=patterns['cutadapt']
    log:
        patterns['cutadapt'] + '.log'
    params:
        extra='-a file:adapters.fa -q 20'
    wrapper:
        wrapper_for('cutadapt')


rule fastqc:
    input: 'samples/{sample}/{sample}{suffix}'
    output:
        html='samples/{sample}/fastqc/{sample}{suffix}_fastqc.html',
        zip='samples/{sample}/fastqc/{sample}{suffix}_fastqc.zip',
    wrapper:
        wrapper_for('fastqc')


rule hisat2:
    input:
        fastq=rules.cutadapt.output.fastq,
        index=[refdict[config['aligner']['tag']]['hisat2']]
    output:
        bam=patterns['bam']
    log:
        patterns['bam'] + '.log'
    wrapper:
        wrapper_for('hisat2/align')


rule fastq_count:
    input:
        fastq='samples/{sample}/{sample}{suffix}.fastq.gz'
    output:
        count='samples/{sample}/{sample}{suffix}.fastq.gz.libsize'
    shell:
        'zcat {input} | echo $((`wc -l`/4)) > {output}'


rule bam_count:
    input:
        bam='samples/{sample}/{sample}{suffix}.bam'
    output:
        count='samples/{sample}/{sample}{suffix}.bam.libsize'
    shell:
        'samtools view -c {input} > {output}'


rule rrna:
    input:
        fastq=rules.cutadapt.output.fastq,
        index=[refdict[config['rrna']['tag']]['bowtie2']]
    output:
        bam=temporary(patterns['rRNA'])
    log:
        patterns['rRNA'] + '.log'
    params: extra='-k 1 --very-fast'
    wrapper:
        wrapper_for('bowtie2/align')


rule featurecounts:
    input:
        annotation=refdict[config['gtf']['tag']]['gtf'],
        bam=rules.hisat2.output
    output:
        counts=patterns['featurecounts']
    log:
        patterns['featurecounts'] + '.log'
    wrapper:
        wrapper_for('featurecounts')


rule libsizes_table:
    input: utils.flatten(targets['libsizes'])
    output: patterns['libsizes_table']
    run:
        def sample(f):
            return os.path.basename(os.path.dirname(f))

        def million(f):
            return float(open(f).read()) / 1e6

        def stage(f):
            return os.path.basename(f).split('.', 1)[1].replace('.gz', '').replace('.count', '')

        df = pd.DataFrame(dict(filename=list(map(str, input))))
        df['sample'] = df.filename.apply(sample)
        df['million'] = df.filename.apply(million)
        df['stage'] = df.filename.apply(stage)
        df = df.set_index('filename')
        df.to_csv(str(output), sep='\t')

rule multiqc:
    input: utils.flatten(targets['fastqc']) + utils.flatten(targets['cutadapt'])
    output: list(set(targets['multiqc']))
    params:
        analysis_directory='samples'
    log: 'multiqc.log'
    wrapper:
        wrapper_for('multiqc')


rule kallisto:
    input:
        index=refdict[config['kallisto']['tag']]['kallisto'],
        fastq=patterns['cutadapt']
    output:
        patterns['kallisto']['h5']
    wrapper:
        wrapper_for('kallisto/quant')


# vim: ft=python
