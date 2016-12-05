import os
import sys
import yaml
import importlib
from snakemake.utils import makedirs
from lcdblib.utils.imports import resolve_name
from lcdblib.utils import utils
from lcdblib.snakemake import aligners, helpers
from common import download_and_postprocess, config_to_dict


HERE = str(srcdir('.'))
def wrapper_for(path):
    return 'file://' + os.path.join(HERE, 'wrappers', 'wrappers', path)


references_dir = os.environ.get('REFERENCES_DIR', config.get('references_dir', None))
if references_dir is None:
    raise ValueError('references dir not set')
config['references_dir'] = references_dir

makedirs([references_dir, os.path.join(references_dir, 'logs')])

references_targets = utils.flatten(config_to_dict(config))

rule all_references:
    input: utils.flatten(config_to_dict(config))


# Downloads the configured URL, applies any configured post-processing, and
# saves the resulting gzipped file to *.fasta.gz or *.gtf.gz.
rule download_and_process:
    output: temporary('{references_dir}/{assembly}/{_type}/{assembly}_{tag}.{_type}.gz')
    run:
        download_and_postprocess(output[0], config, wildcards.assembly, wildcards.tag)


rule unzip:
    input: rules.download_and_process.output
    output: '{references_dir}/{assembly}/{_type}/{assembly}_{tag}.{_type}'
    log: '{references_dir}/logs/{assembly}/{_type}/{assembly}_{tag}.{_type}.log'
    shell: 'gunzip -c {input} > {output}'


rule bowtie2_index:
    output: index=aligners.bowtie2_index_from_prefix('{references_dir}/{assembly}/bowtie2/{assembly}_{tag}')
    input: fasta='{references_dir}/{assembly}/fasta/{assembly}_{tag}.fasta'
    log: '{references_dir}/logs/{assembly}/bowtie2/{assembly}_{tag}.log'
    wrapper: wrapper_for('bowtie2/build')


rule hisat2_index:
    output: index=aligners.hisat2_index_from_prefix('{references_dir}/{assembly}/hisat2/{assembly}_{tag}')
    input: fasta='{references_dir}/{assembly}/fasta/{assembly}_{tag}.fasta'
    log: '{references_dir}/logs/{assembly}/hisat2/{assembly}_{tag}.log'
    wrapper: wrapper_for('hisat2/build')


rule symlink_fasta_to_index_dir:
    input: fasta='{references_dir}/{assembly}/fasta/{assembly}_{tag}.fasta'
    output: '{references_dir}/{assembly}/{index}/{assembly}_{tag}.fasta'
    log: '{references_dir}/logs/{assembly}/{index}/{assembly}_{tag}.fasta.log'
    shell:
        'ln -sf {input} {output}'


rule kallisto_index:
    output: '{references_dir}/{assembly}/kallisto/{assembly}_{tag}.idx'
    input: '{references_dir}/{assembly}/fasta/{assembly}_{tag}.fasta'
    log: '{references_dir}/logs/{assembly}/kallisto/{assembly}_{tag}.log'
    shell:
        '''
        kallisto index -i {output} --make-unique {input} > {log} 2> {log}
        '''

rule conversion_refflat:
    input: '{references_dir}/{assembly}/gtf/{assembly}_{tag}.gtf'
    output: '{references_dir}/{assembly}/gtf/{assembly}_{tag}.refflat'
    log: '{references_dir}/logs/{assembly}/gtf/{assembly}_{tag}.refflat.log'
    conda: 'envs/references_env.yml'
    shell:
        'gtfToGenePred {input} {output}.tmp '
        '''&& awk '{{print $1, $0}}' {output}.tmp > {output} '''
        '&& rm {output}.tmp '


rule chromsizes:
    output: '{references_dir}/{assembly}/fasta/{assembly}_{tag}.chromsizes'
    input: '{references_dir}/{assembly}/fasta/{assembly}_{tag}.fasta'
    log: '{references_dir}/logs/{assembly}/fasta/{assembly}_{tag}.fasta.log'
    conda: 'envs/references_env.yml'
    shell:
        'rm -f {output}.tmp '
        '&& picard CreateSequenceDictionary R={input} O={output}.tmp '
        '&& grep "^@SQ" {output}.tmp '
        '''| awk '{{print $2, $3}}' '''
        '| sed "s/SN://g;s/ LN:/\\t/g" > {output} '
        '&& rm -f {output}.tmp '

# vim: ft=python
