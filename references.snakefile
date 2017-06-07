import os
import sys
import yaml
import importlib
from snakemake.utils import makedirs
from lcdblib.utils.imports import resolve_name
from lcdblib.utils import utils
from lcdblib.snakemake import aligners, helpers
from lib.common import download_and_postprocess, references_dict, get_references_dir

shell.executable('/bin/bash')

localrules: symlink_fasta_to_index_dir, chromsizes

HERE = str(srcdir('.'))

def wrapper_for(path):
    return 'file:' + os.path.join('wrappers', 'wrappers', path)

references_dir = get_references_dir(config)
makedirs([references_dir, os.path.join(references_dir, 'logs')])

refdict, conversion_kwargs = references_dict(config)

rule all_references:
    input: utils.flatten(refdict)


# Downloads the configured URL, applies any configured post-processing, and
# saves the resulting gzipped file to *.fasta.gz or *.gtf.gz.
rule download_and_process:
    output: temporary('{references_dir}/{assembly}/{tag}/{_type}/{assembly}_{tag}.{_type}.gz')
    run:
        download_and_postprocess(output[0], config, wildcards.assembly, wildcards.tag, wildcards._type)


rule unzip:
    input: rules.download_and_process.output
    output: protected('{references_dir}/{assembly}/{tag}/{_type}/{assembly}_{tag}.{_type}')
    log: '{references_dir}/logs/{assembly}/{tag}/{_type}/{assembly}_{tag}.{_type}.log'
    shell: 'gunzip -c {input} > {output}'


rule bowtie2_index:
    output: index=protected(aligners.bowtie2_index_from_prefix('{references_dir}/{assembly}/{tag}/bowtie2/{assembly}_{tag}'))
    input: fasta='{references_dir}/{assembly}/{tag}/fasta/{assembly}_{tag}.fasta'
    log: '{references_dir}/logs/{assembly}/{tag}/bowtie2/{assembly}_{tag}.log'
    wrapper: wrapper_for('bowtie2/build')


rule hisat2_index:
    output: index=protected(aligners.hisat2_index_from_prefix('{references_dir}/{assembly}/{tag}/hisat2/{assembly}_{tag}'))
    input: fasta='{references_dir}/{assembly}/{tag}/fasta/{assembly}_{tag}.fasta'
    log: '{references_dir}/logs/{assembly}/{tag}/hisat2/{assembly}_{tag}.log'
    wrapper: wrapper_for('hisat2/build')


rule symlink_fasta_to_index_dir:
    input: fasta='{references_dir}/{assembly}/{tag}/fasta/{assembly}_{tag}.fasta'
    output: '{references_dir}/{assembly}/{tag}/{index}/{assembly}_{tag}.fasta'
    log: '{references_dir}/logs/{assembly}/{tag}/{index}/{assembly}_{tag}.fasta.log'
    shell:
        'ln -sf {input} {output}'


rule kallisto_index:
    output: protected('{references_dir}/{assembly}/{tag}/kallisto/{assembly}_{tag}.idx')
    input: '{references_dir}/{assembly}/{tag}/fasta/{assembly}_{tag}.fasta'
    log: '{references_dir}/logs/{assembly}/{tag}/kallisto/{assembly}_{tag}.log'
    conda: 'config/envs/references_env.yml'
    shell:
        '''
        kallisto index -i {output} --make-unique {input} > {log} 2> {log}
        '''


rule salmon_index:
    output: protected('{references_dir}/{assembly}/{tag}/salmon/{assembly}_{tag}/hash.bin')
    input:
        fasta='{references_dir}/{assembly}/{tag}/fasta/{assembly}_{tag}.fasta'
    log: '{references_dir}/logs/{assembly}/{tag}/salmon/{assembly}_{tag}.log'
    conda: 'config/envs/references_env.yml'
    wrapper: wrapper_for('salmon/index')


rule conversion_refflat:
    input: '{references_dir}/{assembly}/{tag}/gtf/{assembly}_{tag}.gtf'
    output: protected('{references_dir}/{assembly}/{tag}/gtf/{assembly}_{tag}.refflat')
    log: '{references_dir}/logs/{assembly}/{tag}/gtf/{assembly}_{tag}.refflat.log'
    conda: 'config/envs/references_env.yml'
    shell:
        'gtfToGenePred -ignoreGroupsWithoutExons {input} {output}.tmp '
        '''&& awk '{{print $1"\t"$0}}' {output}.tmp > {output} '''
        '&& rm {output}.tmp '


rule conversion_gffutils:
    input: gtf='{references_dir}/{assembly}/{tag}/gtf/{assembly}_{tag}.gtf'
    output: db=protected('{references_dir}/{assembly}/{tag}/gtf/{assembly}_{tag}.gtf.db')
    log: '{references_dir}/logs/{assembly}/{tag}/gtf/{assembly}_{tag}.gtf.db.log'
    run:
        import gffutils
        kwargs = conversion_kwargs[output[0]]
        db = gffutils.create_db(data=input.gtf, dbfn=output.db, **kwargs)


rule chromsizes:
    output: protected('{references_dir}/{assembly}/{tag}/fasta/{assembly}_{tag}.chromsizes')
    input: '{references_dir}/{assembly}/{tag}/fasta/{assembly}_{tag}.fasta'
    log: '{references_dir}/logs/{assembly}/{tag}/fasta/{assembly}_{tag}.fasta.log'
    conda: 'config/envs/references_env.yml'
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
    input: gtf='{references_dir}/{assembly}/{tag}/gtf/{assembly}_{tag}.gtf'
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
    input: rules.genelist.output
    output:
        protected('{references_dir}/{assembly}/{tag}/gtf/{assembly}_{tag}.{keytype}.csv')
    params:
        prefix='{references_dir}/{assembly}/{tag}/gtf/{assembly}_{tag}',
        ahkey=lambda wildcards, output: conversion_kwargs[output[0]]['ahkey']

    conda: 'config/envs/annotationdbi.yaml'
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
