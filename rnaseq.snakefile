import os
from textwrap import dedent
import yaml
import tempfile
import pandas as pd
from lcdblib.snakemake import helpers, aligners
from lcdblib.utils import utils
from lib import common

# ----------------------------------------------------------------------------
# Note:
#
# In order for automated tests to run, we need to make some settings that are
# not optimal for practical usage. Search this file for the string 
# "# TEST SETTINGS" to find those instances and to learn what changes you might
# want to make.
# ----------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# SETUP
# ----------------------------------------------------------------------------

include: 'references.snakefile'

shell.prefix('set -euo pipefail; export TMPDIR={};'.format(common.tempdir_for_biowulf()))
shell.executable('/bin/bash')
common.get_references_dir(config)

samples, sampletable = common.get_sampletable(config)
refdict, conversion_kwargs = common.references_dict(config)

assembly = config['assembly']

sample_dir = config.get('sample_dir', 'samples')
agg_dir = config.get('aggregation_dir', 'aggregation')

# ----------------------------------------------------------------------------
# PATTERNS
# ----------------------------------------------------------------------------
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
    'rrna_percentages_table': '{agg_dir}/rrna_percentages_table.tsv',
    'rrna_percentages_yaml': '{agg_dir}/rrna_percentages_table_mqc.yaml',
    'rrna': {
        'bam': '{sample_dir}/{sample}/rRNA/{sample}.cutadapt.rrna.bam',
        'libsize': '{sample_dir}/{sample}/rRNA/{sample}.cutadapt.rrna.bam.libsize',
    },
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
    'salmon': '{sample_dir}/{sample}/{sample}.salmon/quant.sf',
    'rseqc': {
        'bam_stat': '{sample_dir}/{sample}/rseqc/{sample}_bam_stat.txt',
    },
    'bigwig': {
        'pos': '{sample_dir}/{sample}/{sample}.cutadapt.bam.pos.bigwig',
        'neg': '{sample_dir}/{sample}/{sample}.cutadapt.bam.neg.bigwig',
    },
    'downstream': {
        'rnaseq': 'downstream/rnaseq.html',
    }
}
fill = dict(sample=samples, sample_dir=sample_dir, agg_dir=agg_dir)
targets = helpers.fill_patterns(patterns, fill)


def wrapper_for(path):
    return 'file:' + os.path.join('wrappers', 'wrappers', path)

# ----------------------------------------------------------------------------
# RULES
# ----------------------------------------------------------------------------

rule targets:
    """
    Final targets to create
    """
    input:
        (
            targets['bam'] +
            utils.flatten(targets['fastqc']) +
            utils.flatten(targets['libsizes']) +
            [targets['fastq_screen']] +
            [targets['libsizes_table']] +
            [targets['rrna_percentages_table']] +
            [targets['multiqc']] +
            utils.flatten(targets['featurecounts']) +
            utils.flatten(targets['rrna']) +
            utils.flatten(targets['markduplicates']) +
            utils.flatten(targets['salmon']) +
            #utils.flatten(targets['dupradar']) +
            utils.flatten(targets['rseqc']) +
            utils.flatten(targets['collectrnaseqmetrics']) +
            utils.flatten(targets['bigwig']) +
            utils.flatten(targets['downstream'])
        )


if 'orig_filename' in sampletable.columns:
    rule symlinks:
        """
        Symlinks files over from original filename
        """
        input: lambda wc: sampletable.set_index(sampletable.columns[0])['orig_filename'].to_dict()[wc.sample]
        output: patterns['fastq']
        run:
            utils.make_relative_symlink(input[0], output[0])

    rule symlink_targets:
        input: targets['fastq']


rule cutadapt:
    """
    Run cutadapt
    """
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
    """
    Run FastQC
    """
    input: '{sample_dir}/{sample}/{sample}{suffix}'
    output:
        html='{sample_dir}/{sample}/fastqc/{sample}{suffix}_fastqc.html',
        zip='{sample_dir}/{sample}/fastqc/{sample}{suffix}_fastqc.zip',
    wrapper:
        wrapper_for('fastqc')


rule hisat2:
    """
    Map reads with HISAT2
    """
    input:
        fastq=rules.cutadapt.output.fastq,
        index=[refdict[assembly][config['aligner']['tag']]['hisat2']]
    output:
        bam=patterns['bam']
    log:
        patterns['bam'] + '.log'
    params:
        samtools_view_extra='-F 0x04'
    threads: 6
    wrapper:
        wrapper_for('hisat2/align')


rule rRNA:
    """
    Map reads with bowtie2 to the rRNA reference
    """
    input:
        fastq=rules.cutadapt.output.fastq,
        index=[refdict[assembly][config['rrna']['tag']]['bowtie2']]
    output:
        bam=patterns['rrna']['bam']
    log:
        patterns['rrna']['bam'] + '.log'
    params:
        samtools_view_extra='-F 0x04'
    threads: 6
    wrapper:
        wrapper_for('bowtie2/align')


rule fastq_count:
    """
    Count reads in a FASTQ file
    """
    input:
        fastq='{sample_dir}/{sample}/{sample}{suffix}.fastq.gz'
    output:
        count='{sample_dir}/{sample}/{sample}{suffix}.fastq.gz.libsize'
    shell:
        'zcat {input} | echo $((`wc -l`/4)) > {output}'


rule bam_count:
    """
    Count reads in a BAM file
    """
    input:
        bam='{sample_dir}/{sample}/{suffix}.bam'
    output:
        count='{sample_dir}/{sample}/{suffix}.bam.libsize'
    shell:
        'samtools view -c {input} > {output}'


rule bam_index:
    """
    Index a BAM
    """
    input:
        bam='{prefix}.bam'
    output:
        bai='{prefix}.bam.bai'
    wrapper: wrapper_for('samtools/index')


rule fastq_screen:
    """
    Run fastq_screen to look for contamination from other genomes
    """
    input:
        fastq=rules.cutadapt.output.fastq,
        dm6=refdict['dmel'][config['aligner']['tag']]['bowtie2'],
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
    """
    Count reads in annotations with featureCounts from the subread package
    """
    input:
        annotation=refdict[assembly][config['gtf']['tag']]['gtf'],
        bam=rules.hisat2.output
    output:
        counts=patterns['featurecounts']
    log:
        patterns['featurecounts'] + '.log'
    wrapper:
        wrapper_for('featurecounts')


rule rrna_libsizes_table:
    """
    Aggregate rRNA counts into a table
    """
    input:
        rrna=targets['rrna']['libsize'],
        fastq=targets['libsizes']['cutadapt']
    output:
        json=patterns['rrna_percentages_yaml'],
        tsv=patterns['rrna_percentages_table']
    run:
        def rrna_sample(f):
            return helpers.extract_wildcards(patterns['rrna']['libsize'], f)['sample']

        def sample(f):
            return helpers.extract_wildcards(patterns['libsizes']['cutadapt'], f)['sample']

        def million(f):
            return float(open(f).read()) / 1e6

        rrna = sorted(input.rrna, key=rrna_sample)
        fastq = sorted(input.fastq, key=sample)
        samples = list(map(rrna_sample, rrna))
        rrna_m = list(map(million, rrna))
        fastq_m = list(map(million, fastq))

        df = pd.DataFrame(dict(
            sample=samples,
            million_reads_rRNA=rrna_m,
            million_reads_fastq=fastq_m,
        ))
        df = df.set_index('sample')
        df['rRNA_percentage'] = df.million_reads_rRNA / df.million_reads_fastq * 100

        df[['million_reads_fastq', 'million_reads_rRNA', 'rRNA_percentage']].to_csv(output.tsv, sep='\t')
        y = {
            'id': 'rrna_percentages_table',
            'section_name': 'rRNA content',
            'description': 'Amount of reads mapping to rRNA sequence',
            'plot_type': 'table',
            'pconfig': {
                'id': 'rrna_percentages_table_table',
                'title': 'rRNA content table',
                'min': 0
            },
            'data': yaml.load(df.transpose().to_json()),
        }
        with open(output.json, 'w') as fout:
            yaml.dump(y, fout, default_flow_style=False)


rule libsizes_table:
    """
    Aggregate fastq and bam counts in to a single table
    """
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
    """
    Aggregate various QC stats and logs into a single HTML report with MultiQC
    """
    input:
        files=(
            utils.flatten(targets['fastqc']) +
            utils.flatten(targets['libsizes_yaml']) +
            utils.flatten(targets['rrna_percentages_yaml']) +
            utils.flatten(targets['cutadapt']) +
            utils.flatten(targets['featurecounts']) +
            utils.flatten(targets['bam']) +
            utils.flatten(targets['markduplicates']) +
            utils.flatten(targets['salmon']) +
            utils.flatten(targets['rseqc']) +
            utils.flatten(targets['fastq_screen']) +
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


rule markduplicates:
    """
    Mark or remove PCR duplicates with Picard MarkDuplicates
    """
    input:
        bam=rules.hisat2.output
    output:
        bam=patterns['markduplicates']['bam'],
        metrics=patterns['markduplicates']['metrics']
    log:
        patterns['markduplicates']['bam'] + '.log'
    params:
        # TEST SETTINGS:
        # You may want to use something larger, like "-Xmx32g" for real-world
        # usage.
        java_args='-Xmx2g'
    wrapper:
        wrapper_for('picard/markduplicates')


rule collectrnaseqmetrics:
    """
    Calculate various RNA-seq QC metrics with Picarc CollectRnaSeqMetrics
    """
    input:
        bam=patterns['bam'],
        refflat=refdict[assembly][config['gtf']['tag']]['refflat']
    output:
        metrics=patterns['collectrnaseqmetrics']['metrics'],
        pdf=patterns['collectrnaseqmetrics']['pdf']
    params:
        # TEST SETTINGS:
        # You may want to use something larger, like "-Xmx32g" for real-world
        # usage.
        java_args='-Xmx32g',
        extra="STRAND=NONE CHART_OUTPUT={}".format(patterns['collectrnaseqmetrics']['pdf'])
    log: patterns['collectrnaseqmetrics']['metrics'] + '.log'
    wrapper: wrapper_for('picard/collectrnaseqmetrics')


rule dupRadar:
    """
    Assess the library complexity with dupRadar
    """
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
    """
    Quantify reads coming from transcripts with Salmon
    """
    input:
        unmatedReads=patterns['cutadapt'],
        index=refdict[assembly][config['salmon']['tag']]['salmon'],
    output: patterns['salmon']
    params: extra="--libType=A"
    log: '{sample_dir}/{sample}/salmon/salmon.quant.log'
    wrapper: wrapper_for('salmon/quant')


rule rseqc_bam_stat:
    """
    Calculate various BAM stats with RSeQC
    """
    input:
        bam=patterns['bam']
    output:
        txt=patterns['rseqc']['bam_stat']
    wrapper: wrapper_for('rseqc/bam_stat')


rule bigwig_neg:
    """
    Create a bigwig for negative-strand reads
    """
    input:
        bam=patterns['bam'],
        bai=patterns['bam'] + '.bai',
    output: patterns['bigwig']['neg']
    threads: 8
    params:
        extra = '--minMappingQuality 20 --ignoreDuplicates --smoothLength 10 --filterRNAstrand forward --normalizeUsingRPKM'
    log:
        patterns['bigwig']['neg'] + '.log'
    wrapper: wrapper_for('deeptools/bamCoverage')


rule bigwig_pos:
    """
    Create a bigwig for postive-strand reads
    """
    input:
        bam=patterns['bam'],
        bai=patterns['bam'] + '.bai',
    output: patterns['bigwig']['pos']
    threads: 8
    params:
        extra = '--minMappingQuality 20 --ignoreDuplicates --smoothLength 10 --filterRNAstrand reverse --normalizeUsingRPKM'
    log:
        patterns['bigwig']['pos'] + '.log'
    wrapper: wrapper_for('deeptools/bamCoverage')


rule rnaseq_rmarkdown:
    """
    Run and render the RMarkdown file that performs differential expression
    """
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
