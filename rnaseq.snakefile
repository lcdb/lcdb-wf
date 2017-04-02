import os
from textwrap import dedent
import yaml
import tempfile
import pandas as pd
from lcdblib.snakemake import helpers, aligners
from lcdblib.utils import utils
from lib import common

TMPDIR = tempfile.gettempdir()
JOBID = os.getenv('SLURM_JOBID')
if JOBID:
    TMPDIR = os.path.join('/lscratch', JOBID)
shell.prefix('set -euo pipefail; export TMPDIR={};'.format(TMPDIR))
shell.executable('/bin/bash')

include: 'references.snakefile'

references_dir = os.environ.get('REFERENCES_DIR', config.get('references_dir', None))
if references_dir is None:
    raise ValueError('No references dir specified')
config['references_dir'] = references_dir

sampletable = pd.read_table(config['sampletable'])
samples = sampletable.ix[:, 0]

assembly = config['assembly']
refdict, conversion_kwargs = common.references_dict(config)

sample_dir = config.get('sample_dir', 'samples')
agg_dir = config.get('aggregation_dir', 'aggregation')

patterns = {
    'fastq':   '{sample_dir}/{sample}/{sample}_R1.fastq.gz',
    'cutadapt': '{sample_dir}/{sample}/{sample}_R1.cutadapt.fastq.gz',
    'bam':     '{sample_dir}/{sample}/{sample}.cutadapt.bam',
    'fastqc': {
        'raw': '{sample_dir}/{sample}/fastqc/{sample}_R1.fastq.gz_fastqc.zip',
        'cutadapt': '{sample_dir}/{sample}/fastqc/{sample}_R1.cutadapt.fastq.gz_fastqc.zip',
        'bam': '{sample_dir}/{sample}/fastqc/{sample}.cutadapt.bam_fastqc.zip',
    },
    'libsizes': {
        'fastq':   '{sample_dir}/{sample}/{sample}_R1.fastq.gz.libsize',
        'cutadapt': '{sample_dir}/{sample}/{sample}_R1.cutadapt.fastq.gz.libsize',
        'bam':     '{sample_dir}/{sample}/{sample}.cutadapt.bam.libsize',
    },
    'fastq_screen': '{sample_dir}/{sample}/{sample}.cutadapt.screen.txt',
    'featurecounts': '{sample_dir}/{sample}/{sample}.cutadapt.bam.featurecounts.txt',
    'libsizes_table': '{agg_dir}/libsizes_table.tsv',
    'libsizes_yaml': '{agg_dir}/libsizes_table_mqc.yaml',
    'multiqc': '{agg_dir}/multiqc.html',
    'markduplicates': {
        'bam': '{sample_dir}/{sample}/{sample}.cutadapt.markdups.bam',
        'metrics': '{sample_dir}/{sample}/{sample}.cutadapt.markdups.bam.metrics',
    },
    'collectrnaseqmetrics': {
        'metrics': '{sample_dir}/{sample}/{sample}.collectrnaseqmetrics.metrics',
        'pdf': '{sample_dir}/{sample}/{sample}.collectrnaseqmetrics.pdf',
    },
    'dupradar': {
        'density_scatter': '{sample_dir}/{sample}/dupradar/{sample}_density_scatter.png',
        'expression_histogram': '{sample_dir}/{sample}/dupradar/{sample}_expression_histogram.png',
        'expression_boxplot': '{sample_dir}/{sample}/dupradar/{sample}_expression_boxplot.png',
        'expression_barplot': '{sample_dir}/{sample}/dupradar/{sample}_expression_barplot.png',
        'multimapping_histogram': '{sample_dir}/{sample}/dupradar/{sample}_multimapping_histogram.png',
        'dataframe': '{sample_dir}/{sample}/dupradar/{sample}_dataframe.tsv',
        'model': '{sample_dir}/{sample}/dupradar/{sample}_model.txt',
        'curve': '{sample_dir}/{sample}/dupradar/{sample}_curve.txt',
    },
    'kallisto': {
        'h5': '{sample_dir}/{sample}/{sample}/kallisto/abundance.h5',
    },
    'salmon': '{sample_dir}/{sample}/{sample}.salmon/quant.sf',
    'rseqc': {
        'bam_stat': '{sample_dir}/{sample}/rseqc/{sample}_bam_stat.txt',
    },
    'downstream': {
        'rnaseq': 'downstream/rnaseq.html',
    }
}
fill = dict(sample=samples, sample_dir=sample_dir, agg_dir=agg_dir)
targets = helpers.fill_patterns(patterns, fill)


def wrapper_for(path):
    return 'file:' + os.path.join('wrappers', 'wrappers', path)

rule targets:
    input:
        (
            targets['bam'] +
            utils.flatten(targets['fastqc']) +
            utils.flatten(targets['libsizes']) +
            [targets['fastq_screen']] +
            [targets['libsizes_table']] +
            [targets['multiqc']] +
            utils.flatten(targets['featurecounts']) +
            utils.flatten(targets['markduplicates']) +
            utils.flatten(targets['dupradar']) +
            utils.flatten(targets['salmon']) +
            utils.flatten(targets['rseqc']) +
            utils.flatten(targets['collectrnaseqmetrics']) +
            utils.flatten(targets['downstream'])
        )


rule cutadapt:
    input:
        fastq=patterns['fastq']
    output:
        fastq=patterns['cutadapt']
    log:
        patterns['cutadapt'] + '.log'
    params:
        extra='-a file:include/adapters.fa -q 20 --minimum-length=25'
    wrapper:
        wrapper_for('cutadapt')


rule fastqc:
    input: '{sample_dir}/{sample}/{sample}{suffix}'
    output:
        html='{sample_dir}/{sample}/fastqc/{sample}{suffix}_fastqc.html',
        zip='{sample_dir}/{sample}/fastqc/{sample}{suffix}_fastqc.zip',
    wrapper:
        wrapper_for('fastqc')


rule hisat2:
    input:
        fastq=rules.cutadapt.output.fastq,
        index=[refdict[assembly][config['aligner']['tag']]['hisat2']]
    output:
        bam=patterns['bam']
    log:
        patterns['bam'] + '.log'
    threads: 6
    wrapper:
        wrapper_for('hisat2/align')


rule fastq_count:
    input:
        fastq='{sample_dir}/{sample}/{sample}{suffix}.fastq.gz'
    output:
        count='{sample_dir}/{sample}/{sample}{suffix}.fastq.gz.libsize'
    shell:
        'zcat {input} | echo $((`wc -l`/4)) > {output}'


rule bam_count:
    input:
        bam='{sample_dir}/{sample}/{sample}{suffix}.bam'
    output:
        count='{sample_dir}/{sample}/{sample}{suffix}.bam.libsize'
    shell:
        'samtools view -c {input} > {output}'


rule fastq_screen:
    input:
        fastq=rules.cutadapt.output.fastq,
        dm6=refdict[assembly][config['aligner']['tag']]['bowtie2'],
        rRNA=refdict[assembly][config['rrna']['tag']]['bowtie2'],
        phix=refdict['phix']['default']['bowtie2']
    output:
        txt=patterns['fastq_screen']
    log:
        patterns['fastq_screen'] + '.log'
    params: subset=100000
    wrapper:
        wrapper_for('fastq_screen')


rule featurecounts:
    input:
        annotation=refdict[assembly][config['gtf']['tag']]['gtf'],
        bam=rules.hisat2.output
    output:
        counts=patterns['featurecounts']
    log:
        patterns['featurecounts'] + '.log'
    wrapper:
        wrapper_for('featurecounts')


rule libsizes_table:
    input:
        utils.flatten(targets['libsizes'])
    output:
        json=patterns['libsizes_yaml'],
        tsv=patterns['libsizes_table']
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
        df = df.pivot('sample', columns='stage', values='million')
        df.to_csv(output.tsv, sep='\t')
        y = {
            'id': 'libsizes_table',
            'section_name': 'Library sizes',
            'description': 'Library sizes at various stages of the pipeline',
            'plot_type': 'table',
            'pconfig': {
                'id': 'libsizes_table_table',
                'title': 'Library size table',
                'min': 0
            },
            'data': yaml.load(df.transpose().to_json()),
        }
        with open(output.json, 'w') as fout:
            yaml.dump(y, fout, default_flow_style=False)


rule multiqc:
    input:
        files=(
            utils.flatten(targets['fastqc']) +
            utils.flatten(targets['libsizes_yaml']) +
            utils.flatten(targets['cutadapt']) +
            utils.flatten(targets['featurecounts']) +
            utils.flatten(targets['bam']) +
            utils.flatten(targets['markduplicates']) +
            utils.flatten(targets['salmon']) +
            utils.flatten(targets['rseqc']) +
            utils.flatten(targets['collectrnaseqmetrics'])
        ),
        config='config/multiqc_config.yaml'
    output: list(set(targets['multiqc']))
    params:
        analysis_directory=" ".join([sample_dir, agg_dir]),
        extra='--config config/multiqc_config.yaml',
    log: list(set(targets['multiqc']))[0] + '.log'
    wrapper:
        wrapper_for('multiqc')


rule kallisto:
    input:
        index=refdict[assembly][config['kallisto']['tag']]['kallisto'],
        fastq=patterns['cutadapt']
    output:
        patterns['kallisto']['h5']
    wrapper:
        wrapper_for('kallisto/quant')


rule markduplicates:
    input:
        bam=rules.hisat2.output
    output:
        bam=patterns['markduplicates']['bam'],
        metrics=patterns['markduplicates']['metrics']
    log:
        patterns['markduplicates']['bam'] + '.log'
    wrapper:
        wrapper_for('picard/markduplicates')


rule collectrnaseqmetrics:
    input:
        bam=patterns['bam'],
        refflat=refdict[assembly][config['gtf']['tag']]['refflat']
    output:
        metrics=patterns['collectrnaseqmetrics']['metrics'],
        pdf=patterns['collectrnaseqmetrics']['pdf']
    params: extra="STRAND=NONE CHART_OUTPUT={}".format(patterns['collectrnaseqmetrics']['pdf'])
    log: patterns['collectrnaseqmetrics']['metrics'] + '.log'
    wrapper: wrapper_for('picard/collectrnaseqmetrics')


rule dupRadar:
    input:
        bam=rules.markduplicates.output.bam,
        annotation=refdict[assembly][config['gtf']['tag']]['gtf'],
    output:
        density_scatter=patterns['dupradar']['density_scatter'],
        expression_histogram=patterns['dupradar']['expression_histogram'],
        expression_boxplot=patterns['dupradar']['expression_boxplot'],
        expression_barplot=patterns['dupradar']['expression_barplot'],
        multimapping_histogram=patterns['dupradar']['multimapping_histogram'],
        dataframe=patterns['dupradar']['dataframe'],
        model=patterns['dupradar']['model'],
        curve=patterns['dupradar']['curve'],
    log: '{sample_dir}/{sample}/dupradar/dupradar.log'
    wrapper:
        wrapper_for('dupradar')


rule salmon:
    input:
        unmatedReads=patterns['cutadapt'],
        index=refdict[assembly][config['salmon']['tag']]['salmon'],
    output: patterns['salmon']
    params: extra="--libType=A"
    log: '{sample_dir}/{sample}/salmon/salmon.quant.log'
    wrapper: wrapper_for('salmon/quant')


rule rseqc_bam_stat:
    input:
        bam=patterns['bam']
    output:
        txt=patterns['rseqc']['bam_stat']
    wrapper: wrapper_for('rseqc/bam_stat')


rule rnaseq_rmarkdown:
    input:
        featurecounts=targets['featurecounts'],
        rmd='downstream/rnaseq.Rmd',
        sampletable=config['sampletable']
    output:
        'downstream/rnaseq.html'
    conda:
        'config/envs/R_rnaseq.yaml'
    shell:
        'Rscript -e '
        '''"rmarkdown::render('{input.rmd}', 'knitrBootstrap::bootstrap_document')"'''

# vim: ft=python
