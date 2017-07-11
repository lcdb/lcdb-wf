import os
from textwrap import dedent
import yaml
import tempfile
import pandas as pd
from lcdblib.snakemake import helpers, aligners
from lcdblib.utils import utils
from lib import common

shell.prefix('set -euo pipefail; export TMPDIR={};'.format(common.tempdir_for_biowulf()))
shell.executable('/bin/bash')

include: 'references.snakefile'

references_dir = os.environ.get('REFERENCES_DIR', config.get('references_dir', None))
if references_dir is None:
    raise ValueError('No references dir specified')
config['references_dir'] = references_dir

samples, sampletable = common.get_sampletable(config)
assembly = config['assembly']
refdict, conversion_kwargs = common.references_dict(config)

sample_dir = config.get('sample_dir', 'samples')
agg_dir = config.get('aggregation_dir', 'aggregation')
fourc_dir = config.get('4c_dir', '4C')

# add to sampletable so we can use zip
sampletable['_sample_dir'] = sample_dir
sampletable['_agg_dir'] = agg_dir
sampletable['_fourc_dir'] = fourc_dir

bait_sequence_lookup = {}
for _, row in sampletable.iterrows():
    key = (row.bait, row.enzyme, str(row.fragLen))
    bait_sequence_lookup[key] = row.primer

patterns = {
    #'bait_stats': '{sample_dir}/{comparison}/{bait}/{bait}_stats.txt',
    'reduced_genome': '{fourc_dir}/{enzyme}_{bait}_flanking_sequences_{fragLen}mer_unique.fa',
    'enzyme': '{fourc_dir}/{enzyme}_site.fa',
    'restriction_sites': '{fourc_dir}/{enzyme}_restriction_sites.bed',
    'flanking_restriction_sites': '{fourc_dir}/{enzyme}_{fragLen}mer_flanking_sites.bed',
    'flanking_restriction_sites_sequences': '{fourc_dir}/{enzyme}_{fragLen}mer_flanking_sites.fa',
    'flanking_restriction_sites_sequences_unique': '{fourc_dir}/{enzyme}_{fragLen}mer_flanking_sites_unique.fa',
    'flanking_restriction_sites_sequences_unique_bed': '{fourc_dir}/{enzyme}_{fragLen}mer_flanking_sites_unique.bed',
    'bowtie2_index': aligners.bowtie2_index_from_prefix('{fourc_dir}/{enzyme}_{fragLen}mer'),
    'fastq': '{sample_dir}/{sample}/{sample}_R1.fastq.gz',
    'bam': '{sample_dir}/{sample}/{enzyme}_{fragLen}mer_trim{fragLen2}/{sample}.bam',
    'unaligned_sam': '{sample_dir}/{sample}/{enzyme}_{fragLen}mer_trim{fragLen2}/{sample}.unaligned.sam',
    'reduced_bam_to_bedgraph': '{sample_dir}/{sample}/{enzyme}_{fragLen}mer_trim{fragLen2}/{sample}.bedgraph',
    'fastqc': {
        'raw': '{sample_dir}/{sample}/fastqc/{sample}_R1.fastq.gz_fastqc.zip',
        #'bam': '{sample_dir}/{sample}/fastqc/{enzyme}_{fragLen}mer_trim{fragLen2}/{sample}.bam_fastqc.zip',
        #'unaligned_sam': '{sample_dir}/{sample}/fastqc/{enzyme}_{fragLen}mer_trim{fragLen2}/{sample}.unaligned.sam_fastq.zip',
    },
    'bait_coords': '{fourc_dir}/{enzyme}_{fragLen}mer_{bait}_bait_coords.bed',
    'multiqc': '{agg_dir}/multiqc.html',
    'remove_bait_and_adjacent_from_bedgraph': '{sample_dir}/{sample}/{enzyme}_{fragLen}mer_trim{fragLen2}/{sample}_cleaned.bedgraph',
    'cleaned_raw_bigwig': '{sample_dir}/{sample}/{enzyme}_{fragLen}mer_trim{fragLen2}/{sample}_cleaned.bigwig',
}



# We typically don't always have the full product of enzyme x sample, so we
# expect all details to be filled out in the sampletable (repeated as needed,
# for example in the case of `enzyme`) and we use the `zip` combination when
# filling in patterns.
targets = helpers.fill_patterns(
    patterns,
    fill=dict(
        sample=sampletable.samplename,
        sample_dir=sampletable._sample_dir,
        agg_dir=sampletable._agg_dir,
        fourc_dir=sampletable._fourc_dir,
        enzyme=sampletable.enzyme,
        bait=sampletable.bait,
        fragLen=sampletable.fragLen,
        fragLen2=sampletable.fragLen2,
    ),
    combination=zip
)

# The next set of patterns are defined according to the "4c" section of the
# config.
#
patterns_4cker = {
    '4cker':
        '4cker-output/{comparison}/sentinel.txt',
    'cis_bedgraphs':
        '4cker-output/{comparison}/cis_k{cis_k}/{sample}_cis_norm_counts.bedGraph',
    'cis_bigwigs':
        '4cker-output/{comparison}/cis_k{cis_k}/{sample}_cis_norm_counts.bigwig',
    'nearbait_bedgraphs':
        '4cker-output/{comparison}/nearbait_k{nearbait_k}/{sample}_nearbait_norm_counts.bedGraph',
    'nearbait_bigwigs':
        '4cker-output/{comparison}/nearbait_k{nearbait_k}/{sample}_nearbait_norm_counts.bigwig',
    'cis_adaptive_bed':
        '4cker-output/{comparison}/cis_k{cis_k}/{bait}_cis_adaptive_windows.bed',
    'cis_adaptive_bigbed':
        '4cker-output/{comparison}/cis_k{cis_k}/{bait}_cis_k{cis_k}_adaptive_windows.bigbed',
    'nearbait_adaptive_bed':
        '4cker-output/{comparison}/nearbait_k{nearbait_k}/{bait}_nearbait_adaptive_windows.bed',
    'nearbait_adaptive_bigbed':
        '4cker-output/{comparison}/nearbait_k{nearbait_k}/{bait}_nearbait_adaptive_windows.bigbed',
    'nearbait_sample_interacting_bigbed':
        '4cker-output/{comparison}/nearbait_k{nearbait_k}/{sample}_nearbait_{inter}.bigbed',
    'nearbait_treatment_interacting_bigbed':
        '4cker-output/{comparison}/nearbait_k{nearbait_k}/{bait}_{treatment}_nearbait_highinter.bigbed',
    'cis_sample_interacting_bigbed':
        '4cker-output/{comparison}/cis_k{cis_k}/{sample}_cis_{inter}.bigbed',
    'cis_treatment_interacting_bigbed':
        '4cker-output/{comparison}/cis_k{cis_k}/{bait}_{treatment}_cis_highinter.bigbed',
    'cis_colorized': '4cker-output/{comparison}/cis_k{cis_k}/{bait}_cis_colorized_differential.bigbed',
    'nearbait_colorized': '4cker-output/{comparison}/nearbait_k{nearbait_k}/{bait}_nearbait_colorized_differential.bigbed',
}

fills = []
for comparison, vals in config['4c']['comparisons'].items():
    bait = vals['bait']
    d = {
        'comparison': comparison,
        'bait': bait,
        'cis_k': config['4c']['baits'][bait]['cis_k'],
        'nearbait_k': config['4c']['baits'][bait]['nearbait_k'],
    }
    for treatment in ['control', 'treatment']:
        for sample in vals[treatment]:
            for inter in ['noninter', 'lowinter', 'highinter']:
                d['inter'] = inter
                d['samplename'] = sample
                d['treatment'] = treatment
                fills.append(d.copy())

fills = pd.DataFrame(fills)

targets_4cker = helpers.fill_patterns(
    patterns_4cker,
    fill=dict(
        comparison=fills.comparison,
        bait=fills.bait,
        cis_k=fills.cis_k,
        sample=fills.samplename,
        nearbait_k=fills.nearbait_k,
        treatment=fills.treatment,
        inter=fills.inter,
    ),
    combination=zip
)

patterns.update(patterns_4cker)
targets.update(targets_4cker)


def wrapper_for(path):
    return 'file:' + os.path.join('wrappers', 'wrappers', path)

rule targets:
    input:
        (
            utils.flatten(targets['cleaned_raw_bigwig']) +
            utils.flatten(targets['bait_coords']) +
            utils.flatten(targets['4cker']) +
            utils.flatten(targets['multiqc'])

        )

rule trackhub:
    input:
        (
            utils.flatten(targets['cis_bigwigs']) +
            utils.flatten(targets['cis_adaptive_bigbed']) +
            utils.flatten(targets['nearbait_bigwigs']) +
            utils.flatten(targets['nearbait_adaptive_bigbed']) +
            utils.flatten(targets['cis_treatment_interacting_bigbed']) +
            utils.flatten(targets['cis_sample_interacting_bigbed']) +
            utils.flatten(targets['nearbait_treatment_interacting_bigbed']) +
            utils.flatten(targets['nearbait_sample_interacting_bigbed']) +
            utils.flatten(targets['cis_colorized']) +
            utils.flatten(targets['nearbait_colorized'])
        )

rule all_fastqc:
    input: utils.flatten(targets['fastqc'])

rule fastqc:
    input: '{sample_dir}/{sample}/{sample}{suffix}'
    output:
        html='{sample_dir}/{sample}/fastqc/{sample}{suffix}_fastqc.html',
        zip='{sample_dir}/{sample}/fastqc/{sample}{suffix}_fastqc.zip',
    wrapper:
        wrapper_for('fastqc')

# ----------------------------------------------------------------------------
# Create a FASTA file for the enzyme's restriction site
rule restriction_seq:
    output: patterns['enzyme']
    run:
        from Bio import Restriction
        with open(output[0], 'w') as fout:
            fout.write(
                '> {0}\n{1}'.format(
                    wildcards.enzyme,
                    getattr(Restriction, wildcards.enzyme).site
                )
            )


# ----------------------------------------------------------------------------
# BED file of restriction sites in the genome
rule restriction_sites:
    input:
        sequence=refdict[assembly][config['4c']['tag']]['fasta'],
        oligos=patterns['enzyme']
    output:
        patterns['restriction_sites']
    params:
        enzyme='{enzyme}'
    wrapper: wrapper_for('oligomatch')


# ----------------------------------------------------------------------------
# Flanking regions on either side of restriction sites. Does not include the
# site itself.
rule flanking_restriction_sites:
    input:
        sites=rules.restriction_sites.output,
        chromsizes=refdict[assembly][config['4c']['tag']]['chromsizes'],
    output: patterns['flanking_restriction_sites']
    conda: 'config/envs/4c.yaml'
    shell:
        'bedtools flank '
        '-i {input.sites} '
        '-g {input.chromsizes} '
        '-b {wildcards.fragLen} '
        '| sort -k1,1 -k2,2n -k3,3n '
        '> {output}'


# ----------------------------------------------------------------------------
# Sequences of flanking regions
rule flanking_restriction_sites_sequences:
    input:
        genome=refdict[assembly][config['4c']['tag']]['fasta'],
        sites=rules.flanking_restriction_sites.output
    output: patterns['flanking_restriction_sites_sequences']
    conda: 'config/envs/4c.yaml'
    shell:
        'bedtools getfasta '
        '-fi {input.genome} '
        '-bed {input.sites} '
        '-fo {output} '

# See https://github.com/rr1859/R.4Cker/issues/26
# ----------------------------------------------------------------------------
# Remove multimapping sequences
# rule flanking_restriction_sites_sequences_unique:
#     input: rules.flanking_restriction_sites_sequences.output
#     output: patterns['flanking_restriction_sites_sequences_unique']
#     shell:
#         "grep -v '^>' {input} "          # get just the sequences
#         "| sort | uniq -i -u "           # case-insensitive; duplicated seqs are omitted
#         "| grep -xF -f - -B 1 {input} "  # exact matching seq from input, one per line, plus the coords in the line above
#         "| grep -v '^--' > {output}"     # remove `--` output by -B1

rule flanking_restriction_sites_sequences_unique:
    input: rules.flanking_restriction_sites_sequences.output
    output: patterns['flanking_restriction_sites_sequences_unique']
    shell:
        "paste - - < {input} "
        "| awk '{{a[toupper($2)]++; b[$2] = $1}}; END "
        """{{for(n in a) if (a[n] == 1) print b[n], n}}' """
        "| sed 's/-\|:/\\t/g' "
        "| sort -k1,1 -k2,2n -k3,3n "
        """| awk '{{print $1":"$2"-"$3"\\n"$4}}' """
        "> {output} "


# ----------------------------------------------------------------------------
# Convert the header names of uniqued sequences into BED
# E.g.,
#
#   >chr2L:1-100
#
# becomes:
#
#   chr2L    1    100
#
rule flanking_restriction_sites_sequences_unique_bed:
    input: rules.flanking_restriction_sites_sequences_unique.output
    output: patterns['flanking_restriction_sites_sequences_unique_bed']
    shell:
        "grep '^>' {input} "
        "| sed 's/>//g' "
        "| sed 's/:\|-/\\t/g' "
        "> {output}"

# ----------------------------------------------------------------------------
# Build Bowtie2 index on the reduced genome
rule bowtie2_index_reduced_genome:
    output:
        index=patterns['bowtie2_index']
    input:
        fasta=rules.flanking_restriction_sites_sequences_unique.output
    log: patterns['bowtie2_index'][0] + '.log'
    wrapper: wrapper_for('bowtie2/build')


# ----------------------------------------------------------------------------
# Map fastqs to reduced genome

class Default(dict):
    def __missing__(self, key):
        return '{' + key + '}'


rule bowtie2_map_to_reduced_genome:
    input:
        index=(
            rules.bowtie2_index_reduced_genome.output[0]
            .format_map(Default(fourc_dir=fourc_dir))
        ),
        fastq=patterns['fastq'],
    output:
        bam=patterns['bam'],
        unaligned_sam=patterns['unaligned_sam']
    threads: 8
    params:
        bowtie2_extra=lambda wildcards, output:
            ('-N 0 -5 {0} --un {1}'.format(wildcards.fragLen2, output.unaligned_sam))
    log: patterns['bam'] + '.bowtie2.log'
    wrapper: wrapper_for('bowtie2/align')


# ----------------------------------------------------------------------------
# Converts a BAM file with sequence names of the form "chrom:start-stop" into
# a bedGraph file
rule reduced_bam_to_bedgraph:
    input: rules.bowtie2_map_to_reduced_genome.output.bam
    output: patterns['reduced_bam_to_bedgraph']
    conda: 'config/envs/4c.yaml'
    shell:
        'samtools view {input} '
        "| awk '{{print $3}}' "             # sequence column of BAM
        '| grep -v "^*" '                   # ignore unmapped
        '| sort | uniq -c '                 # reads in each fragment
        "| sed -e 's/\:\| \|-/\\t/g' "      # convert coord to tab-delimited BED3
        '''| awk '{{OFS="\\t"; print $2, $3, $4, $1}}' '''  # bedGraph format
        '> {output}'


# ----------------------------------------------------------------------------
# Identifies the coordinates of the bait fragment and the immediately-adjacent
# fragments as a BED file.
rule identify_bait_coords:
    input:
        fasta=rules.flanking_restriction_sites_sequences_unique.output,
        bed=rules.flanking_restriction_sites_sequences_unique_bed.output,
    output: patterns['bait_coords']
    run:
        from Bio import SeqIO
        from Bio import Seq

        sequence = bait_sequence_lookup[(wildcards.bait, wildcards.enzyme, wildcards.fragLen)]

        res = []
        primer_forward = sequence.lower()
        primer_reverse = str(Seq.Seq(sequence).reverse_complement()).lower()
        for i, rec in enumerate(SeqIO.parse(open(input.fasta[0]), 'fasta')):
            rec_seq = rec.seq.lower()
            if (primer_forward in rec_seq) or (primer_reverse in rec_seq):
                coord = rec.name
                chrom, startstop = coord.split(':')
                start, stop = startstop.split('-')
                res.append((i, (chrom, start, stop)))
        if len(res) == 0:
            raise ValueError("Primer not found")
        if len(res) > 1:
            raise ValueError("More than one primer found")

        # fragment number and coordinates of the bait
        pos, coords = res[0]

        # now we find
        with open(output[0], 'w') as fout:
            for i, line in enumerate(open(input.bed[0])):

                # if we're at the position of the bait, one before, or one
                # after then print them
                if (pos - 1) <= i <= (pos + 1):

                    # We found this line by its order in the file; double-check
                    # that its coords match that of the found bait from above.
                    if i == pos:
                        assert tuple(line.strip().split('\t')[:3]) == coords
                    fout.write(line)
                continue


# ----------------------------------------------------------------------------
# Remove the bait and the immediately adjacent sequences from the bedgraph
def _bait_coords_for_sample(wildcards):
    """
    Extracts information from the sampletable to fill in the correct bait in
    the filename, which is otherwise missing from the output.
    """
    mapping = list(
        sampletable[sampletable.samplename == wildcards.sample]
        .transpose()
        .to_dict().values()
    )[0]
    return patterns['bait_coords'].format(fourc_dir=fourc_dir, **mapping)

rule remove_bait_and_adjacent_from_bedgraph:
    input:
        bedgraph=rules.reduced_bam_to_bedgraph.output,
        coords=_bait_coords_for_sample
    output: patterns['remove_bait_and_adjacent_from_bedgraph']
    shell:
        'bedtools intersect -v -a {input.bedgraph} -b {input.coords} | grep -v -E "Un|random" > {output}'


# ----------------------------------------------------------------------------
# Aggregate fastqc
rule multiqc:
    input:
        files=(
            utils.flatten(targets['fastqc'])
            + utils.flatten(targets['bam'])
        ),
        config='config/multiqc_config.yaml'
    output: list(set(targets['multiqc']))
    params:
        analysis_directory=" ".join([sample_dir]),
        extra='--config config/multiqc_config.yaml',
    log: list(set(targets['multiqc']))[0] + '.log'
    wrapper:
        wrapper_for('multiqc')



rule R4cker:
    input:
        bedgraphs=utils.flatten(targets['remove_bait_and_adjacent_from_bedgraph']),
        config='config/4c-config.yml',
        rscript='downstream/4c.R',
    output: touch('4cker-output/{comparison}/sentinel.txt')
    log: '4cker-output/{comparison}/log'
    shell:
        'source activate 4c-wf '
        '&& Rscript --default-packages=boot,cluster,foreign,lattice,MASS,Matrix,nlme,base,compiler,datasets,graphics,grDevices,grid,methods,parallel,splines,stats,stats4,tcltk,tools,utils '
        '{input.rscript} --config {input.config} --comparison {wildcards.comparison} &> {log}'

rule bedgraph_to_bigwig:
    input:
        chromsizes=refdict[assembly][config['4c']['tag']]['chromsizes'],
        bedgraph='4cker-output/{comparison}/{kind}/{prefix}.bedGraph'
    output:
        bigwig='4cker-output/{comparison}/{kind}/{prefix}.bigwig'
    run:
        comparison = wildcards.comparison
        bait = config['4c']['comparisons'][comparison]['bait']
        chrom = config['4c']['baits'][bait]['chrom']
        windowsize = 1000
        shell(
            'python lib/fourc/bedgraph_to_bigwig.py '
            '--chromosome {chrom} '
            '--chromsizes {input.chromsizes} '
            '--windowsize {windowsize} '
            '--bedgraph {input.bedgraph} '
            '--output {input.bedgraph}.windowed'
        )
        shell(
            'bedGraphToBigWig '
            '{input.bedgraph}.windowed '
            '{input.chromsizes} '
            '{output.bigwig} '
        )

rule raw_bedgraph_to_bigwig:
    input:
        chromsizes=refdict[assembly][config['4c']['tag']]['chromsizes'],
        bedgraph='{sample_dir}/{sample}/{enzyme}_{fragLen}mer_trim{fragLen2}/{sample}_cleaned.bedgraph',
    output:
        '{sample_dir}/{sample}/{enzyme}_{fragLen}mer_trim{fragLen2}/{sample}_cleaned.bigwig',
    shell:
        'bedtools sort -i {input.bedgraph} | '
        'bedtools merge -i - -c 4 -o sum > {input.bedgraph}.sorted '
        '&& bedGraphToBigWig '
        '{input.bedgraph}.sorted '
        '{input.chromsizes} '
        '{output} '
        '&& rm {input.bedgraph}.sorted'


def _colorized_input(wildcards):
    comparison = wildcards.comparison
    bait = config['4c']['comparisons'][comparison]['bait']
    kind = wildcards.kind
    k = config['4c']['baits'][bait][kind + '_k']
    control_samples = config['4c']['comparisons'][comparison]['control']
    treatment_samples = config['4c']['comparisons'][comparison]['treatment']
    return utils.flatten({
        'control': expand(
            '4cker-output/{comparison}/{kind}_k{k}/{sample}_{kind}_norm_counts.bedGraph',
            comparison=comparison, kind=kind, k=k, sample=control_samples, bait=bait),
        'treatment': expand(
            '4cker-output/{comparison}/{kind}_k{k}/{sample}_{kind}_norm_counts.bedGraph',
            comparison=comparison, kind=kind, k=k, sample=treatment_samples, bait=bait),
        'diff': expand(
            '4cker-output/{comparison}/{kind}_diff_k{k}/{bait}_control_treatment_{kind}_pval0.05_diff.bed',
            comparison=comparison, kind=kind, k=k, sample=treatment_samples, bait=bait),
    })


rule colorized:
    input: _colorized_input
    output: '4cker-output/{comparison}/{kind}_k{k}/{bait}_{kind}_colorized_differential.bed'
    run:
        comparison = wildcards.comparison
        bait = wildcards.bait
        kind = wildcards.kind
        k = wildcards.k
        control_samples = config['4c']['comparisons'][comparison]['control']
        treatment_samples = config['4c']['comparisons'][comparison]['treatment']
        def u(x):
            return sorted(list(set(x)))
        files = {
            'control': ' '.join(u(expand(
                '4cker-output/{comparison}/{kind}_k{k}/{sample}_{kind}_norm_counts.bedGraph',
                comparison=comparison, kind=kind, k=k, sample=control_samples, bait=bait))),
            'treatment': ' '.join(u(expand(
                '4cker-output/{comparison}/{kind}_k{k}/{sample}_{kind}_norm_counts.bedGraph',
                comparison=comparison, kind=kind, k=k, sample=treatment_samples, bait=bait))),
            'diff': ' '.join(u(expand(
                '4cker-output/{comparison}/{kind}_diff_k{k}/{bait}_control_treatment_{kind}_pval0.05_diff.bed',
                comparison=comparison, kind=kind, k=k, sample=treatment_samples, bait=bait))),
        }
        shell(
            'python lib/fourc/find-up-dn.py '
            '--bed {files[diff]} '
            '--output {output} '
            '--control {files[control]} '
            '--treatment {files[treatment]} '
        )

rule bed_to_bigbed:
    input:
        bed='{prefix}.bed',
        chromsizes=refdict[assembly][config['4c']['tag']]['chromsizes'],
    output: '{prefix}.bigbed'
    run:
        if 'colorized' in output[0]:
            autosql = dedent("""
            table lfc
            "Browser extensible data (9 fields) with float lfc."
                (
                string chrom;      "Chromosome (or contig, scaffold, etc.)"
                uint   chromStart; "Start position in chromosome"
                uint   chromEnd;   "End position in chromosome"
                string name;       "Name of item"
                uint score;        "score"
                char[1] strand;    "+ or -"
                uint thickStart;   "Start of where display should be thick (start codon)"
                uint thickEnd;     "End of where display should be thick (stop codon)"
                uint reserved;     "Used as itemRgb as of 2004-11-22"
                string log2foldchange; "Log2 fold change"
                )
            """)
            with open(output[0] + '.as', 'w') as fout:
                fout.write(autosql)
            as_arg = '-type=bed9+ -as={}.as '.format(output[0])
        else:
            as_arg = ''

        shell(
            """
            if [ -s {input.bed} ]; then
                GENOME=$(mktemp);
                awk '{{OFS="\\t"; print $1, "0", $2}}' {input.chromsizes} > $GENOME
                bedtools intersect -a <(sed "s/ /\\t/g" {input.bed}) -b $GENOME > ${{GENOME}}.tmp
                bedToBigBed {as_arg} ${{GENOME}}.tmp {input.chromsizes} {output}
                rm $GENOME ${{GENOME}}.tmp
            else
                touch {output}
            fi
            """
        )

# vim: ft=python
