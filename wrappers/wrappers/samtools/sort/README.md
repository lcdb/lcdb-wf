# Wrapper for samtools sort.

## Example:

```
rule samtools_sort:
    input:
        "mapped/{sample}.bam"
    output:
        "mapped/{sample}.sorted.bam"
    params:
        extra: "-m 4G"
    threads: 8
    wrapper:
        "file://path/to/wrapper"
```
