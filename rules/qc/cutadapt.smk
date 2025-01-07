rule cutadapt:
    input:
        fastq=expand(patterns["fastq"], n=n, allow_missing=True)
    output:
        fastq=expand(patterns["cutadapt"], n=n, allow_missing=True)
    log:
        'data/rnaseq_samples/{sample}/{sample}_cutadapt.fastq.gz.log'
    threads: 6
    resources:
        mem_mb=gb(2),
        runtime=autobump(hours=2)
    params:
        extra=(
            "--nextseq-trim 20 "
            "--overlap 6 "
            "--minimum-length 25 "
            "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA "
        ) + "-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT " if is_paired else ""
    run:
        if is_paired:
            shell(
                "cutadapt "
                "-o {output[0]} "
                "-p {output[1]} "
                "-j {threads} "
                "{params.extra} "
                "{input.fastq[0]} "
                "{input.fastq[1]} "
                "&> {log}"
            )
        else:
            shell(
                "cutadapt "
                "-o {output[0]} "
                "-j {threads} "
                "{params.extra} "
                "{input.fastq[0]} "
                "&> {log}"
            )
