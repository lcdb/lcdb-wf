
rule sample_strand_check:
    input:
        fastq=fill_r1_r2(c.sampletable, c.patterns['fastq']),
        index=rules.bowtie2_index.output,
        bed12=rules.conversion_bed12.output,
    output:
        strandedness='strand_check/{sample}/{sample}.strandedness',
        bam=temporary('strand_check/{sample}/{sample}.strandedness.bam'),
        bai=temporary('strand_check/{sample}/{sample}.strandedness.bam.bai'),
        fastqs=temporary(expand('strand_check/{sample}/{sample}_R{n}.strandedness.fastq', sample=SAMPLES, n=n)),
    log:
        'strand_check/{sample}/{sample}.strandedness.log'
    threads: 6
    resources:
        mem_mb=gb(8),
        runtime=autobump(hours=2)
    run:
        prefix = aligners.prefix_from_bowtie2_index(input.index)
        nreads = int(config['strand_check_reads']) * 4
        if c.is_paired:
            assert len(input.fastq) == 2
            assert len(output.fastqs) == 2
            shell('set +o pipefail; zcat {input.fastq[0]} | head -n {nreads} > {output.fastqs[0]}')
            shell('set +o pipefail; zcat {input.fastq[0]} | head -n {nreads} > {output.fastqs[1]}')
            fastqs = f'-1 {output.fastqs[0]} -2 {output.fastqs[1]} '
        else:
            assert len(input.fastq) == 1
            assert len(output.fastqs) == 1
            shell('set +o pipefail; zcat {input.fastq[0]} | head -n {nreads} > {output.fastqs[0]}')
            fastqs = f'-U {output.fastqs[0]} '
        shell(
            "bowtie2 "
            "-x {prefix} "
            "{fastqs} "
            '--no-unal '
            "--threads {threads} 2> {log} "
            "| samtools view -Sb - "
            "| samtools sort - -o {output.bam} "
        )
        shell("samtools index {output.bam}")
        shell(
            'infer_experiment.py -r {input.bed12} -i {output.bam} > {output} 2> {log}'
        )

rule strand_check:
    input:
        expand('strand_check/{sample}/{sample}.strandedness', sample=SAMPLES)
    output:
        html='strand_check/strandedness.html',
        filelist=temporary('strand_check/filelist')
    log:
        'strand_check/strandedness.log'
    resources:
        mem_mb=gb(1),
        runtime=autobump(hours=2)
    run:
        with open(output.filelist, 'w') as fout:
            for i in  input:
                fout.write(i + '\n')
        shell(
            'multiqc '
            '--force '
            '--module rseqc '
            '--file-list {output.filelist} '
            '--filename {output.html} &> {log}'
        )

# vim: ft=snakemake
