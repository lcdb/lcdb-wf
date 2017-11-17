import os
import sys
sys.path.insert(0, srcdir('../..'))
import yaml
import importlib
import tempfile
from snakemake.utils import makedirs
from lcdblib.utils.imports import resolve_name
from lcdblib.utils import utils
from lcdblib.snakemake import aligners, helpers
from lib import common

tempfile.tempdir = common.tempdir_for_biowulf()

shell.prefix('set -euo pipefail; export TMPDIR={};'.format(common.tempdir_for_biowulf()))
shell.executable('/bin/bash')
references_dir = common.get_references_dir(config)
refdict, conversion_kwargs = common.references_dict(config)

makedirs([references_dir, os.path.join(references_dir, 'logs')])


def wrapper_for(path):
    return 'file:' + os.path.join('../..','wrappers', 'wrappers', path)


localrules: symlink_fasta_to_index_dir, chromsizes

rule all_references:
    input: utils.flatten(refdict)


rule download_and_process:
    """Downloads the configured URL, applies any configured post-processing, and
    saves the resulting gzipped file to *.fasta.gz or *.gtf.gz.
    """
    output:
        temporary('{references_dir}/{assembly}/{tag}/{_type}/{assembly}_{tag}.{_type}.gz')
    run:
        common.download_and_postprocess(output[0], config, wildcards.assembly, wildcards.tag, wildcards._type)


rule unzip:
    """Generic rule to unzip files as needed, for example when building
    indexes.
    """
    input:
        rules.download_and_process.output
    output:
        protected('{references_dir}/{assembly}/{tag}/{_type}/{assembly}_{tag}.{_type}')
    log:
        '{references_dir}/logs/{assembly}/{tag}/{_type}/{assembly}_{tag}.{_type}.log'
    shell: 'gunzip -c {input} > {output}'


rule bowtie2_index:
    """
    Build bowtie2 index
    """
    input:
        '{references_dir}/{assembly}/{tag}/fasta/{assembly}_{tag}.fasta'
    output:
        protected(aligners.bowtie2_index_from_prefix('{references_dir}/{assembly}/{tag}/bowtie2/{assembly}_{tag}'))
    log:
        '{references_dir}/logs/{assembly}/{tag}/bowtie2/{assembly}_{tag}.log'
    run:
        prefix = aligners.prefix_from_bowtie2_index(output)
        shell(
            'bowtie2-build '
            '{input} '
            '{prefix} '
            '&> {log}')


rule hisat2_index:
    """
    Build HISAT2 index
    """
    input:
        '{references_dir}/{assembly}/{tag}/fasta/{assembly}_{tag}.fasta'
    output:
        protected(aligners.hisat2_index_from_prefix('{references_dir}/{assembly}/{tag}/hisat2/{assembly}_{tag}'))
    log:
        '{references_dir}/logs/{assembly}/{tag}/hisat2/{assembly}_{tag}.log'
    run:
        prefix = aligners.prefix_from_hisat2_index(output)
        shell(
            'hisat2-build '
            '{input} '
            '{prefix} '
            '&> {log}')


rule symlink_fasta_to_index_dir:
    """Aligners often want the reference fasta in the same dir as the index, so
    this makes the appropriate symlink
    """
    input:
        fasta='{references_dir}/{assembly}/{tag}/fasta/{assembly}_{tag}.fasta'
    output:
        '{references_dir}/{assembly}/{tag}/{index}/{assembly}_{tag}.fasta'
    log:
        '{references_dir}/logs/{assembly}/{tag}/{index}/{assembly}_{tag}.fasta.log'
    shell:
        'ln -sf {input} {output}'


rule salmon_index:
    "Build salmon index"
    output:
        protected('{references_dir}/{assembly}/{tag}/salmon/{assembly}_{tag}/hash.bin')
    input:
        fasta='{references_dir}/{assembly}/{tag}/fasta/{assembly}_{tag}.fasta'
    log:
        '{references_dir}/logs/{assembly}/{tag}/salmon/{assembly}_{tag}.log'
    params:
        outdir='{references_dir}/{assembly}/{tag}/salmon/{assembly}_{tag}'
    shell:
        'salmon index '
        '--transcripts {input.fasta} '
        '--index {params.outdir} '
        '&> {log}'


rule conversion_refflat:
    """Converts a GTF into refFlat format
    """
    input:
        '{references_dir}/{assembly}/{tag}/gtf/{assembly}_{tag}.gtf'
    output:
        protected('{references_dir}/{assembly}/{tag}/gtf/{assembly}_{tag}.refflat')
    log:
        '{references_dir}/logs/{assembly}/{tag}/gtf/{assembly}_{tag}.refflat.log'
    shell:
        'gtfToGenePred -ignoreGroupsWithoutExons {input} {output}.tmp '
        '''&& awk '{{print $1"\t"$0}}' {output}.tmp > {output} '''
        '&& rm {output}.tmp '


rule conversion_gffutils:
    """Converts a GTF into a gffutils sqlite3 database
    """
    input:
        gtf='{references_dir}/{assembly}/{tag}/gtf/{assembly}_{tag}.gtf'
    output:
        db=protected('{references_dir}/{assembly}/{tag}/gtf/{assembly}_{tag}.gtf.db')
    log:
        '{references_dir}/logs/{assembly}/{tag}/gtf/{assembly}_{tag}.gtf.db.log'
    run:
        import gffutils
        kwargs = conversion_kwargs[output[0]]
        fd, tmpdb = tempfile.mkstemp(suffix='.db', prefix='gffutils_')
        db = gffutils.create_db(data=input.gtf, dbfn=tmpdb, **kwargs)
        shell('mv {tmpdb} {output.db}')


rule chromsizes:
    """Creates a chromsizes table from fasta
    """
    input:
        '{references_dir}/{assembly}/{tag}/fasta/{assembly}_{tag}.fasta'
    output:
        protected('{references_dir}/{assembly}/{tag}/fasta/{assembly}_{tag}.chromsizes')
    log:
        '{references_dir}/logs/{assembly}/{tag}/fasta/{assembly}_{tag}.fasta.log'
    shell:
        'export LC_COLLATE=C; '
        'rm -f {output}.tmp '
        '&& picard CreateSequenceDictionary R={input} O={output}.tmp '
        '&& grep "^@SQ" {output}.tmp '
        '''| awk '{{print $2, $3}}' '''
        '| sed "s/SN://g;s/ LN:/\\t/g" '
        '| sort -k1,1 > {output} '
        '&& rm -f {output}.tmp '


rule genelist:
    """Creates a list of unique gene names in the GTF
    """
    input:
        gtf='{references_dir}/{assembly}/{tag}/gtf/{assembly}_{tag}.gtf'
    output:
        protected('{references_dir}/{assembly}/{tag}/gtf/{assembly}_{tag}.genelist')
    run:
        attribute = conversion_kwargs[output[0]]['gene_id']
        import gffutils
        genes = set()
        for feature in gffutils.DataIterator(input.gtf):
            genes.update(feature.attributes[attribute])
        with open(output[0], 'w') as fout:
            for feature in sorted(list(set(genes))):
                fout.write(feature + '\n')


rule annotations:
    """Creates TSVs mappings between gene names (created in the `genelist` rule
    above) and all of the columns in the configured AnnotationHub accession
    (https://bioconductor.org/packages/release/bioc/html/AnnotationHub.html)
    """
    input:
        rules.genelist.output
    output:
        protected('{references_dir}/{assembly}/{tag}/gtf/{assembly}_{tag}.{keytype}.csv')
    params:
        prefix='{references_dir}/{assembly}/{tag}/gtf/{assembly}_{tag}',
        ahkey=lambda wildcards, output: conversion_kwargs[output[0]]['ahkey']
    shell:
        '''Rscript -e "'''
        "library(AnnotationHub); "
        "ah <- AnnotationHub(); "
        "db <- ah[['{params.ahkey}']]; "
        "gene.names <- read.table('{input}', stringsAsFactors=FALSE)[['V1']];"
        "for (col in columns(db)){{"
        "f <- select(db, keys=gene.names, keytype='{wildcards.keytype}', columns=col);"
        "write.csv(f, file=paste0('{params.prefix}', '.', col, '.csv'), row.names=FALSE);"
        '''}}"'''


# vim: ft=python
