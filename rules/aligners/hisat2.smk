rule hisat2:
    input:
        fastq=expand(patterns["cutadapt"], n=n, allow_missing=True),
        index=rules.hisat2_index.output,
    output:
        bam=temporary(patterns['bam'])
    log:
        patterns['bam'] + '.log'
    threads: 6
    resources:
        mem_mb=gb(32),
        runtime=autobump(hours=8)
    params:
        extra=""
    run:
        prefix = os.path.commonprefix(input.index).rstrip(".")
        sam = output.bam.replace('.bam', '.sam')

        if is_paired:
            assert len(input.fastq) == 2
            fastqs = '-1 {0} -2 {1} '.format(*input.fastq)
        else:
            assert len(input.fastq) == 1
            fastqs = '-U {0} '.format(input.fastq)

        shell(
            "hisat2 "
            "-x {prefix} "
            "{fastqs} "
            '--no-unal '
            "--threads {threads} "
            "-S {sam} "
            "> {log} 2>&1"
        )

        shell(
            "samtools view -Sb {sam} "
            "| samtools sort - -o {output.bam} -O BAM "
            "&& rm {sam}"
        )
