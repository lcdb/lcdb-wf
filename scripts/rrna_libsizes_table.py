"""
Prepares a TSV and JSON file for multiqc to pick up and display as a sortable
table
"""
import os
import re
import pandas as pd
import yaml
from snakemake.io import regex_from_filepattern


def rrna_sample(f):
    m = re.compile(
        regex_from_filepattern(
            snakemake.params.rrna_pattern,
        )
    ).match(f)
    if m:
        return m.groupdict()["sample"]


def sample(f):
    m = re.compile(
        regex_from_filepattern(
            snakemake.params.fastq_pattern,
        )
    ).match(f)
    if m:
        return m.groupdict()["sample"]


def million(f):
    return float(open(f).read()) / 1e6


rrna = sorted(snakemake.input.rrna, key=rrna_sample)
fastq = sorted(snakemake.input.fastq, key=sample)
samples = list(map(rrna_sample, rrna))
rrna_m = list(map(million, rrna))
fastq_m = list(map(million, fastq))

df = pd.DataFrame(
    dict(
        sample=samples,
        million_reads_rRNA=rrna_m,
        million_reads_fastq=fastq_m,
    )
)
df = df.set_index("sample")
df["rRNA_percentage"] = df.million_reads_rRNA / df.million_reads_fastq * 100

df[["million_reads_fastq", "million_reads_rRNA", "rRNA_percentage"]].to_csv(
    snakemake.output.tsv, sep="\t"
)
y = {
    "id": "rrna_percentages_table",
    "section_name": "rRNA content",
    "description": "Amount of reads mapping to rRNA sequence",
    "plot_type": "table",
    "pconfig": {
        "id": "rrna_percentages_table_table",
        "title": "rRNA content table",
        "min": 0,
    },
    "data": yaml.load(df.transpose().to_json(), Loader=yaml.FullLoader),
}
with open(snakemake.output.json, "w") as fout:
    yaml.dump(y, fout, default_flow_style=False)
