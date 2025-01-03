import os
import sys
import pandas

HERE = str(Path(workflow.snakefile).parent)
sys.path.insert(0, HERE + "/../..")
from lib.utils import autobump, gb, hours
from lib import utils

REFERENCES = config.get('reference_dir', '../../references')

def default_postprocess(origfn, newfn):
    shell("mv {origfn} {newfn}")

rule fasta:
    output:
        temporary(REFERENCES + '/genome.fa.gz')
    log:
        REFERENCES + "/logs/genome.fa.gz.log"
    run:
        utils.download_and_postprocess(
            urls=config['fasta']['url'],
            postprocess=config['fasta'].get('postprocess', None),
            outfile=output[0],
            log=log
        )


rule gtf:
    output:
        temporary(REFERENCES + '/annotation.gtf.gz')
    log:
        REFERENCES + "/logs/annotation.gtf.gz.log"
    run:
        utils.download_and_postprocess(
            urls=config['gtf']['url'],
            postprocess=config['gtf'].get('postprocess', None),
            outfile=output[0],
            log=log
        )


rule rrna:
    output:
        temporary(REFERENCES + '/rrna.fa.gz')
    log:
        REFERENCES + "/logs/rrna.fa.gz.log"
    run:
        utils.download_and_postprocess(
            urls=config['rrna']['url'],
            postprocess=config['rrna'].get('postprocess', None),
            outfile=output[0],
            log=log
        )


rule unzip:
    input:
        REFERENCES + '/{prefix}.gz'
    output:
        REFERENCES + '/{prefix}'
    shell: 'gunzip -c {input} > {output}'


rule bowtie2_index:
    input:
        REFERENCES + '/{label}.fa',
    output:
        multiext(
            REFERENCES + '/bowtie2/{label}',
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
            ".fa",
        ),
    log:
        REFERENCES + '/logs/bowtie2_{label}.log'
    resources:
        runtime=autobump(hours=8),
        mem_mb=autobump(gb=32),
        disk_mb=autobump(gb=50)
    threads:
        8
    run:
        index = os.path.commonprefix(output).rstrip(".")
        shell(
            "bowtie2-build"
            " --threads {threads}"
            " {input}"
            " {index}"
            " &> {log}"
        )
        utils.make_relative_symlink(input[0], output[-1])


rule star_index:
    input:
        fasta=REFERENCES + '/genome.fa',
        gtf=REFERENCES + '/annotation.gtf',
    output:
        REFERENCES + '/star/Genome'
    log:
        REFERENCES + '/logs/star.log'
    threads:
        8
    resources:
        runtime=autobump(hours=8),
        mem_mb=gb(64)
    run:
        genomedir = os.path.dirname(output[0])
        shell('rm -r {genomedir}')
        shell('mkdir -p {genomedir}')
        shell(
            'STAR '
            '--runMode genomeGenerate '
            '--runThreadN {threads} '
            '--genomeDir {genomedir} '
            '--genomeFastaFiles {input.fasta} '

            # NOTE: GTF is optional
            '--sjdbGTFfile {input.gtf} '

            # NOTE: STAR docs say that 100 should work well.
            '--sjdbOverhang 100 '

            # NOTE: for small genomes, may need to scale this down to
            # min(14, log2(GenomeLength) / 2 - 1)
            # --genomeSAindexNbases 14
            '&> {log}'
        )
        # STAR writes a hard-coded Log.out file to the current working
        # directory. So put that on the end of the log file for the rule and
        # then clean up.
        shell('cat {genomedir}/Log.out >> {log} && rm {genomedir}/Log.out')
        shell("ln -s {input.fasta} {genomedir}")

rule hisat2_index:
    input:
        REFERENCES + '/genome.fa',
    output:
        multiext(
            REFERENCES + '/hisat2/genome',
            ".1.ht2",
            ".2.ht2",
            ".3.ht2",
            ".4.ht2",
            ".5.ht2",
            ".6.ht2",
            ".7.ht2",
            ".8.ht2",
            ".fa",
        )
    log:
        REFERENCES + '/logs/hisat2.log'
    resources:
        runtime=autobump(hours=8),
        mem_mb=autobump(gb=32),
        disk_mb=autobump(gb=50)
    threads:
        8
    run:
        index = os.path.commonprefix(output).rstrip(".")
        shell(
            "hisat2-build"
            " --threads {threads}"
            " {input}"
            " {index}"
            " &> {log}"
        )
        shell("ln -s {input} {output[-1]}")



rule transcriptome_fasta:
    input:
        fasta=REFERENCES + '/genome.fa',
        gtf=REFERENCES + '/annotation.gtf',
    output:
        REFERENCES + '/transcriptome.fa' 
    resources:
        runtime=hours(1)
    shell:
        'gffread {input.gtf} -w {output} -g {input.fasta}'


rule salmon_index:
    input:
        REFERENCES + '/transcriptome.fa'
    output:
        REFERENCES + '/salmon/versionInfo.json'
    log:
        REFERENCES + '/logs/salmon.log'
    params:
        outdir=REFERENCES + '/salmon'
    resources:
        mem_mb=gb(32),
        runtime=hours(2)
    run:
        outdir = os.path.dirname(output[0])
        shell(
            'salmon index '
            '--transcripts {input} '
            '--index {outdir} '
            '&> {log}'
        )


rule kallisto_index:
    output:
        REFERENCES + '/kallisto/transcripts.idx',
    input:
        REFERENCES + '/genome.fa'
    log:
        REFERENCES + '/logs/kallisto.log'
    resources:
        runtime=hours(2),
        mem_mb=gb(32),
    shell:
        'kallisto index '
        '--index {output} '
        '{input} '
        '&> {log}'


rule conversion_refflat:
    input:
        REFERENCES + '/annotation.gtf'
    output:
        REFERENCES + '/annotation.refflat'
    log:
        REFERENCES + '/logs/annotation.refflat.log'
    resources:
        runtime=hours(2),
        mem_mb=gb(2)
    shell:
        'gtfToGenePred -ignoreGroupsWithoutExons {input} {output}.tmp '
        '''&& awk '{{print $1"\t"$0}}' {output}.tmp > {output} '''
        '&& rm {output}.tmp '


rule conversion_bed12:
    input:
        REFERENCES + '/annotation.gtf'
    output:
        REFERENCES + '/annotation.bed12'
    resources:
        runtime=hours(2),
        mem_mb=gb(2)
    shell:
        'gtfToGenePred -ignoreGroupsWithoutExons {input} {output}.tmp '
        '&& genePredToBed {output}.tmp {output} '
        '&& rm {output}.tmp'


rule chromsizes:
    input:
        REFERENCES + '/genome.fa'
    output:
        REFERENCES + '/genome.chromsizes'
    log:
        REFERENCES + '/logs/genome.chromsizes.log'
    params:
        # NOTE: Be careful with the memory here; make sure you have enough
        # and/or it matches the resources you're requesting
        java_args='-Xmx20g'
        # java_args='-Xmx2g'  # [TEST SETTINGS -1]
    resources:
        mem_mb=gb(24),
        runtime=hours(2)
    shell:
        'export LC_COLLATE=C; '
        'rm -f {output}.tmp '
        '&& picard '
        '{params.java_args} '
        'CreateSequenceDictionary R={input} O={output}.tmp &> {log} '
        '&& grep "^@SQ" {output}.tmp '
        '''| awk '{{print $2, $3}}' '''
        '| sed "s/SN://g;s/ LN:/\\t/g" '
        '| sort -k1,1 > {output} '
        '&& rm -f {output}.tmp '


rule mappings:
    """
    Creates gzipped TSV mapping between attributes in the GTF.
    """
    input:
        gtf=REFERENCES + '/annotation.gtf'
    output:
        REFERENCES + '/annotation.mapping.tsv.gz'
    params:
        include_featuretypes=lambda wildcards, output: conversion_kwargs[output[0]].get('include_featuretypes', [])
    resources:
        runtime=hours(2),
        mem_mb=gb(2)
    run:
        import gffutils

        # Will want to change the setting back to what it was originally when
        # we're done
        orig_setting = gffutils.constants.always_return_list
        gffutils.constants.always_return_list = False

        include_featuretypes = params.include_featuretypes

        res = []
        for f in gffutils.DataIterator(input[0]):

            ft = f.featuretype

            if include_featuretypes and (ft not in include_featuretypes):
                continue

            d = dict(f.attributes)
            d['__featuretype__'] = ft
            res.append(d)

        df = pandas.DataFrame(res)

        # Depending on how many attributes there were and the
        # include_featuretypes settings, this may take a while.
        df = df.drop_duplicates()

        df.to_csv(output[0], sep='\t', index=False, compression='gzip')

        # Restore original setting
        gffutils.constants.always_return_list = orig_setting
