# Wrapper for samtools view

## Example:

```
rule samtools_view:
    input:
        bam="mapped/{sample}.bam"
    output:
        bam="mapped/{sample}.sorted.bam"
    params:
        extra: "-F 0x04"
    wrapper:
        "file://path/to/wrapper"
```
