import pandas as pd
from lcdblib.snakemake import helpers, aligners
from lcdblib.utils import utils

shell.prefix('set -euo pipefail; ')

include: 'references/references.snakefile'

sampletable = pd.read_table(config['sampletable'])
samples = sampletable.ix[:, 0]


patterns = {
    'fastq':   'samples/{sample}/{sample}_R1.fastq.gz',
    'trimmed': 'samples/{sample}/{sample}_R1.trim.fastq.gz',
    'bam':     'samples/{sample}/{sample}.trim.bam',
    'fastqc': {
        'raw': 'samples/{sample}/fastqc/{sample}_R1.fastq.gz.html',
        'trim': 'samples/{sample}/fastqc/{sample}_R1.trim.fastq.gz.html',
        'bam': 'samples/{sample}/fastqc/{sample}.trim.bam.html',
    },
    'counts': {
        'fastq':   'samples/{sample}/{sample}_R1.fastq.gz.count',
        'trimmed': 'samples/{sample}/{sample}_R1.trim.fastq.gz.count',
        'bam':     'samples/{sample}/{sample}.trim.bam.count',
        'rRNA':     'samples/{sample}/{sample}.trim.rRNA.bam.count',
    },
    'rRNA': 'samples/{sample}/{sample}.trim.rRNA.bam',
    'featurecounts': 'samples/{sample}/{sample}.trim.bam.featurecounts.txt',
    'counts_table': 'counts_table.tsv',
    'multiqc': 'multiqc.html'
}
fill = dict(sample=samples, count=['.count', ''])
targets = helpers.fill_patterns(patterns, fill)


HERE = str(srcdir(''))
def wrapper_for(path):
    return 'file://' + os.path.join(HERE, 'wrappers', 'wrappers', path)


rule targets:
    input: targets['bam'] + utils.flatten(targets['fastqc']) + utils.flatten(targets['counts']) + [targets['counts_table']] + [targets['multiqc']]


rule cutadapt:
    input:
        fastq=patterns['fastq']
    output:
        fastq=patterns['trimmed']
    log:
        patterns['fastq'] + '.log'
    params:
        extra='-a file:adapters.fa -q 20'
    wrapper:
        wrapper_for('cutadapt')


rule fastqc:
    input: 'samples/{sample}/{sample}{suffix}'
    output:
        html='samples/{sample}/fastqc/{sample}{suffix}.html',
        zip='samples/{sample}/fastqc/{sample}{suffix}.zip',
    wrapper:
        wrapper_for('fastqc')


rule hisat2:
    input:
        fastq=rules.cutadapt.output.fastq,
        index=aligners.hisat2_index_from_prefix(config['aligner_prefix'])
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
        count='samples/{sample}/{sample}{suffix}.fastq.gz.count'
    shell:
        'zcat {input} | echo $((`wc -l`/4)) > {output}'


rule bam_count:
    input:
        bam='samples/{sample}/{sample}{suffix}.bam'
    output:
        count='samples/{sample}/{sample}{suffix}.bam.count'
    shell:
        'samtools view -c {input} > {output}'


rule rrna:
    input:
        fastq=rules.cutadapt.output.fastq,
        index=aligners.bowtie2_index_from_prefix(config['rrna_prefix'])
    output:
        bam=temporary(patterns['rRNA'])
    log:
        patterns['rRNA'] + '.log'
    params: extra='-k 1 --very-fast'
    wrapper:
        wrapper_for('bowtie2/align')


rule featurecounts:
    input:
        annotation=config['gtf'],
        bam=rules.hisat2.output
    output:
        counts=patterns['featurecounts']
    log:
        patterns['featurecounts'] + '.log'
    wrapper:
        wrapper_for('featurecounts')


rule counts_table:
    input: utils.flatten(targets['counts'])
    output: patterns['counts_table']
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
    input: utils.flatten(targets['fastqc']) + [i + '.log' for i in utils.flatten(targets['cutadapt'])]
    output: targets['multiqc']
    log: 'multiqc.log'
    wrapper:
        wrapper_for('multiqc')

# vim: ft=python
