import os
import sys
import yaml
import importlib
from lcdblib.utils.imports import resolve_name
from lcdblib.snakemake import aligners, helpers
from common import download_and_postprocess


def wrapper_for(path):
    return os.path.join('file://', str(srcdir('.')), '..', 'wrappers', 'wrappers', path)


references_dir = config['references_dir']
if not os.path.exists(references_dir):
    os.makedirs(references_dir)


# Map "type" in config to file extensions
ext_mapping = {
    'gtf': '.gtf',
    'fasta': '.fa.gz',
}

# Map "indexes" value to a pattern specific to each index.
index_mapping = {
    'bowtie2': ('{references_dir}/{assembly}/bowtie2/{assembly}{tag}.{n}.bt2', dict(n=[1, 2, 3, 4])),
    'hisat2': ('{references_dir}/{assembly}/hisat2/{assembly}{tag}.{n}.ht2', dict(n=range(1, 9))),
    'kallisto': ('{references_dir}/{assembly}/kallisto/{assembly}{tag}.idx', dict()),
}

references_targets = []

for block in config['references']:

    # build up the local vars for filling in the pattern
    tag = '_' + block.get('tag', 'default')
    ext = ext_mapping[block['type']]
    assembly = block['assembly']

    references_targets.append('{references_dir}/{assembly}/{assembly}{tag}{ext}'.format(**locals()))

    # Add the indexes too, if they were asked for
    if block['type'] == 'fasta':
        indexes = block.get('indexes', [])
        for index in indexes:
            pattern, kwargs = index_mapping[index]
            kwargs = kwargs.copy()
            references_targets.extend(
                expand(
                    pattern, assembly=assembly, tag=tag,
                    references_dir=references_dir, **kwargs
                )
            )


rule all_references:
    input: references_targets

# The downloading rules all have the same general form and support arbitrary
# post-processing functions to be specified in the config file.

rule download_fasta:
    output: '{references_dir}/{assembly}/{assembly}_{tag}.fa.gz'
    run:
        download_and_postprocess(output[0], config, wildcards.assembly, wildcards.tag)


rule download_rrna:
    output: '{references_dir}/{assembly}/rRNA/{assembly}_{tag}.fa.gz'
    run:
        download_and_postprocess(output[0], config, wildcards.assembly, wildcards.tag)


rule unzip_fasta:
    input: rules.download_fasta.output
    output: temporary('{references_dir}/{assembly}/{assembly}_{tag}.fa')
    shell: 'gunzip -c {input} > {output}'


rule download_gtf:
    output: '{references_dir}/{assembly}/{assembly}_{tag}.gtf.gz'
    run:
        download_and_postprocess(output[0], config, wildcards.assembly, wildcards.tag)


rule unzip_gtf:
    input: rules.download_gtf.output
    output: '{references_dir}/{assembly}/{assembly}_{tag}.gtf'
    shell: 'gunzip -c {input} > {output}'

rule bowtie_index:
    output: aligners.bowtie2_index_from_prefix('{references_dir}/{assembly}/bowtie2/{assembly}_{tag}')
    input: rules.unzip_fasta.output
    log: '{references_dir}/{assembly}/bowtie2/{assembly}{tag}.log'
    shell:
        '''
        bowtie2-build {input} {references_dir}/{wildcards.assembly}/bowtie2/{wildcards.assembly}_{wildcards.tag} > {log} 2> {log}
        '''


rule hisat2_index:
    input:
        fasta=rules.unzip_fasta.output
    output:
        index=aligners.hisat2_index_from_prefix('{references_dir}/{assembly}/hisat2/{assembly}_{tag}')
    log:
        '{references_dir}/{assembly}/hisat2/{assembly}{tag}.log'
    wrapper: wrapper_for('hisat2/build')


rule kallisto_index:
    output: '{references_dir}/{assembly}/kallisto/{assembly}_{tag}.idx'
    input: '{references_dir}/{assembly}/{assembly}{tag}.fa.gz'
    log: '{references_dir}/{assembly}/kallisto/{assembly}{tag}.log'
    shell:
        '''
        kallisto index -i {output} --make-unique {input} > {log} 2> {log}
        '''
# vim: ft=python
