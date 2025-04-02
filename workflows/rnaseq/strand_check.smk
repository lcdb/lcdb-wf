rule sample_strand_check:
    input:
        fastq=expand("data/rnaseq_samples/{{sample}}/{{sample}}_R{n}.fastq.gz", n=n),
        index=expand(rules.bowtie2_index.output, label="genome"),
        bed12=rules.conversion_bed12.output,
    output:
        strandedness="strand_check/{sample}/{sample}.strandedness",
        bam=temporary("strand_check/{sample}/{sample}.strandedness.bam"),
        bai=temporary("strand_check/{sample}/{sample}.strandedness.bam.bai"),
        fastqs=temporary(
            expand(
                "strand_check/{sample}/{sample}_R{n}.strandedness.fastq",
                n=n,
                allow_missing=True,
            )
        ),
    log:
        "strand_check/{sample}/{sample}.strandedness.log",
    threads: 6
    resources:
        mem="8g",
        runtime="2h",
    run:
        prefix = os.path.commonprefix(input.index).rstrip(".")
        nreads = int(1e5 * 4)
        if is_paired:
            shell(
                "set +o pipefail; zcat {input.fastq[0]} | head -n {nreads} > {output.fastqs[0]}"
            )
            shell(
                "set +o pipefail; zcat {input.fastq[0]} | head -n {nreads} > {output.fastqs[1]}"
            )
            fastqs = f"-1 {output.fastqs[0]} -2 {output.fastqs[1]} "
        else:
            shell(
                "set +o pipefail; zcat {input.fastq[0]} | head -n {nreads} > {output.fastqs[0]}"
            )
            fastqs = f"-U {output.fastqs[0]} "
        shell(
            "bowtie2 "
            "-x {prefix} "
            "{fastqs} "
            "--no-unal "
            "--threads {threads} 2> {log} "
            "| samtools view -Sb - "
            "| samtools sort - -o {output.bam} "
        )
        shell("samtools index {output.bam}")
        shell(
            "infer_experiment.py -r {input.bed12} -i {output.bam} > {output} 2> {log}"
        )


rule strand_check:
    input:
        expand("strand_check/{sample}/{sample}.strandedness", sample=SAMPLES),
    output:
        html="strand_check/strandedness.html",
        filelist=temporary("strand_check/filelist"),
    log:
        "strand_check/strandedness.log",
    resources:
        mem="1g",
        runtime="2h",
    run:
        with open(output.filelist, "w") as fout:
            for i in input:
                fout.write(i + "\n")
        shell(
            "multiqc "
            "--force "
            "--module rseqc "
            "--file-list {output.filelist} "
            "--filename {output.html} &> {log}"
        )
