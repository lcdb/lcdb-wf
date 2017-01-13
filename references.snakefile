import os
import sys
import yaml
import importlib
from snakemake.utils import makedirs
from lcdblib.utils.imports import resolve_name
from lcdblib.utils import utils
from lcdblib.snakemake import aligners, helpers
from common import download_and_postprocess, references_dict, get_references_dir


HERE = str(srcdir('.'))
def wrapper_for(path):
    return 'file://' + os.path.join('wrappers', 'wrappers', path)


references_dir = get_references_dir(config)
makedirs([references_dir, os.path.join(references_dir, 'logs')])


rule all_references:
    input: utils.flatten(references_dict(config))


# Downloads the configured URL, applies any configured post-processing, and
# saves the resulting gzipped file to *.fasta.gz or *.gtf.gz.
rule download_and_process:
    output: temporary('{references_dir}/{assembly}/{tag}/{_type}/{assembly}_{tag}.{_type}.gz')
    run:
        download_and_postprocess(output[0], config, wildcards.assembly, wildcards.tag, wildcards._type)


rule unzip:
    input: rules.download_and_process.output
    output: '{references_dir}/{assembly}/{tag}/{_type}/{assembly}_{tag}.{_type}'
    log: '{references_dir}/logs/{assembly}/{tag}/{_type}/{assembly}_{tag}.{_type}.log'
    shell: 'gunzip -c {input} > {output}'


rule bowtie2_index:
    output: index=aligners.bowtie2_index_from_prefix('{references_dir}/{assembly}/{tag}/bowtie2/{assembly}_{tag}')
    input: fasta='{references_dir}/{assembly}/{tag}/fasta/{assembly}_{tag}.fasta'
    log: '{references_dir}/logs/{assembly}/{tag}/bowtie2/{assembly}_{tag}.log'
    wrapper: wrapper_for('bowtie2/build')


rule hisat2_index:
    output: index=aligners.hisat2_index_from_prefix('{references_dir}/{assembly}/{tag}/hisat2/{assembly}_{tag}')
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
    output: '{references_dir}/{assembly}/{tag}/kallisto/{assembly}_{tag}.idx'
    input: '{references_dir}/{assembly}/{tag}/fasta/{assembly}_{tag}.fasta'
    log: '{references_dir}/logs/{assembly}/{tag}/kallisto/{assembly}_{tag}.log'
    conda: 'envs/references_env.yml'
    shell:
        '''
        kallisto index -i {output} --make-unique {input} > {log} 2> {log}
        '''

rule conversion_refflat:
    input: '{references_dir}/{assembly}/{tag}/gtf/{assembly}_{tag}.gtf'
    output: '{references_dir}/{assembly}/{tag}/gtf/{assembly}_{tag}.refflat'
    log: '{references_dir}/logs/{assembly}/{tag}/gtf/{assembly}_{tag}.refflat.log'
    conda: 'envs/references_env.yml'
    shell:
        'gtfToGenePred -ignoreGroupsWithoutExons {input} {output}.tmp '
        '''&& awk '{{print $1, $0}}' {output}.tmp > {output} '''
        '&& rm {output}.tmp '


rule chromsizes:
    output: '{references_dir}/{assembly}/{tag}/fasta/{assembly}_{tag}.chromsizes'
    input: '{references_dir}/{assembly}/{tag}/fasta/{assembly}_{tag}.fasta'
    log: '{references_dir}/logs/{assembly}/{tag}/fasta/{assembly}_{tag}.fasta.log'
    conda: 'envs/references_env.yml'
    shell:
        'rm -f {output}.tmp '
        '&& picard CreateSequenceDictionary R={input} O={output}.tmp '
        '&& grep "^@SQ" {output}.tmp '
        '''| awk '{{print $2, $3}}' '''
        '| sed "s/SN://g;s/ LN:/\\t/g" > {output} '
        '&& rm -f {output}.tmp '

# vim: ft=python
