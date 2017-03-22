# Wrapper for `kallisto quant`


[Kallisto home](https://pachterlab.github.io/kallisto/)

[Kallisto manual](https://pachterlab.github.io/kallisto/manual)

## Examples

Single-end reads example, with recommended `extra` args:

```python
rule kallisto_quant:
    input:
        fastq="{sample}.fastq.gz",
        index="references/dm6/kallisto.idx"
    output:
        h5="quant/kallisto/{sample}/abundance.h5",
        tsv="quant/kallisto/{sample}/abundance.tsv",
        json="quant/kallisto/{sample}/run_info.json"
    log: "quant/kallisto/{sample}/kallisto.log"
    params: extra="--fragment-length=200 --sd=20 --single"
    wrapper:
        "/path/to/wrapper/location"
```

Paired-end, simplified example:

```python
rule kallisto_quant:
    input:
        fastq=["{sample}_R1.fastq.gz", "{sample}_R2.fastq.gz"]
        index="references/dm6/kallisto.idx"
    output:
        h5="quant/kallisto/{sample}/abundance.h5",
        tsv="quant/kallisto/{sample}/abundance.tsv",
        json="quant/kallisto/{sample}/run_info.json"
```

## Input

* `index`: kallisto index file

* `fastq`: either a pair of PE fastqs or a single fastq. If single, see note in *Params* section.

## Output

- `h5`: HDF5 quantification file
- `tsv`: TSV quantification file
- `json`: JSON metadata

`kallisto quant` outputs files hard-coded as:

```
outdir/
├── abundance.h5
├── abundance.tsv
└── run_info.json
```

Since the downstream `sleuth` tool expects these hard-coded names, we do not
move them. It is recommended that you use the sample name in the directory path
to separate them (see examples).

The `--output-dir` argument to `kallisto quant` is determined from the output
files, and the wrapper checks to make sure they are all in the same directory.

## Threads

Threads are passed as the `--threads` argument to `kallisto quant`.

## Params

* `extra`: passed verbatim to `kallisto quant`

*NOTE:* if you are using single-end reads, then kallisto will complain if you
don't set the `--fragment-length` and `--sd` params. Add them to extra. Based
on [this
comment](https://groups.google.com/forum/#!searchin/kallisto-sleuth-users/sd/kallisto-sleuth-users/VPJfzL502bw/e2JDq7ezBgAJ)
from the author, `--fragment-length=200 --sd=20` should be a good starting
point.

*NOTE:* if you plan on using `sleuth` for differential expression, you'll want
to include the `--bootstrap-samples` argument.

