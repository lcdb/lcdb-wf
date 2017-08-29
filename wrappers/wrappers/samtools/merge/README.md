# Wrapper for samtools merge

## Example:

```
rule samtools_merge:
    input:
        "mapped/{sample}_1.bam",
        "mapped/{sample}_2.bam"
    output:
        "mapped/{sample}.merged.bam"
    params:
        extra: "-m 4G"
    threads: 8
    wrapper:
        "file://path/to/wrapper"
```
