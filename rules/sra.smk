
sampletable['orig_filename'] = expand(
    'original_data/sra_samples/{sample}/{sample}_R{n}.fastq.gz', sample=SAMPLES, n=1)

if is_paired:
    sampletable['orig_filename_R2'] = expand(
        'original_data/sra_samples/{sample}/{sample}_R{n}.fastq.gz', sample=SAMPLES, n=2)

rule fastq_dump:
    output:
        fastq=expand('original_data/sra_samples/{sample}/{sample}_R{n}.fastq.gz', n=n)
    log:
        'original_data/sra_samples/{sample}/{sample}.fastq.gz.log'
    params:
        is_paired=is_paired,
        sampletable=_st,
        # extra="-X 100000",  # [TEST SETTINGS]
    resources:
        mem_mb=gb(1),
        disk_mb=autobump(gb=1),
        runtime=autobump(hours=2)
    run:
        _st = sampletable.set_index(sampletable.columns[0])
        srr = _st.loc[wildcards.sample, "Run"]
        extra = params.get("extra", "")
        if is_paired:
            shell("fastq-dump {srr} --gzip --split-files {extra} &> {log}")
            shell("mv {srr}_1.fastq.gz {output[0]}")
            shell("mv {srr}_2.fastq.gz {output[1]}")
        else:
            shell("fastq-dump {srr} -Z {extra} 2> {log} | gzip -c > {output[0]}.tmp")
            shell("mv {output[0]}.tmp {output[0]}")

# vim: ft=snakemake
