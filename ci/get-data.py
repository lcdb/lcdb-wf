#!/usr/bin/env python
import argparse
import os

from snakemake.shell import shell
from snakemake.utils import makedirs

BRANCH = "master"
URL = "https://github.com/lcdb/lcdb-test-data/blob/{0}/data/{{}}?raw=true".format(
    BRANCH
)

TOPLEVEL = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _download_file(fn, dest=None, verbose=False):
    url = URL.format(fn)
    if dest is None:
        dest = fn
    dest = os.path.join(TOPLEVEL, dest)
    makedirs(os.path.dirname(dest))
    if not verbose:
        q = "-q"
    else:
        q = ""
    shell(f"wget {q} -O- {url} > {dest}")
    if verbose:
        print(f"Saved {dest}")
    return dest


ap = argparse.ArgumentParser()
ap.add_argument("-v", "--verbose", action="store_true", help="Be verbose when downloading")
args = ap.parse_args()

_download_file(
    "rnaseq_samples/sample1/sample1.small_R1.fastq.gz",
    "workflows/rnaseq/data/example_data/rnaseq_sample1.fq.gz",
    args.verbose,
)
_download_file(
    "rnaseq_samples/sample2/sample2.small_R1.fastq.gz",
    "workflows/rnaseq/data/example_data/rnaseq_sample2.fq.gz",
    args.verbose,
)
_download_file(
    "rnaseq_samples/sample3/sample3.small_R1.fastq.gz",
    "workflows/rnaseq/data/example_data/rnaseq_sample3.fq.gz",
    args.verbose,
)
_download_file(
    "rnaseq_samples/sample4/sample4.small_R1.fastq.gz",
    "workflows/rnaseq/data/example_data/rnaseq_sample4.fq.gz",
    args.verbose,
)

_download_file(
    "rnaseq_samples/sample1/sample1.small_R1.fastq.gz",
    "workflows/rnaseq/data/example_data/rnaseq_sample1PE_1.fq.gz",
    args.verbose,
)
_download_file(
    "rnaseq_samples/sample1/sample1.small_R2.fastq.gz",
    "workflows/rnaseq/data/example_data/rnaseq_sample1PE_2.fq.gz",
    args.verbose,
)
_download_file(
    "rnaseq_samples/sample2/sample2.small_R1.fastq.gz",
    "workflows/rnaseq/data/example_data/rnaseq_sample2PE_1.fq.gz",
    args.verbose,
)
_download_file(
    "rnaseq_samples/sample2/sample2.small_R2.fastq.gz",
    "workflows/rnaseq/data/example_data/rnaseq_sample2PE_2.fq.gz",
    args.verbose,
)

_download_file(
    "chipseq_samples/input_1/input_1.tiny_R1.fastq.gz",
    "workflows/chipseq/data/example_data/chipseq_input1.fq.gz",
    args.verbose,
)
_download_file(
    "chipseq_samples/ip_1/ip_1.tiny_R1.fastq.gz",
    "workflows/chipseq/data/example_data/chipseq_ip1.fq.gz",
    args.verbose,
)
_download_file(
    "chipseq_samples/input_2/input_2.tiny_R1.fastq.gz",
    "workflows/chipseq/data/example_data/chipseq_input2.fq.gz",
    args.verbose,
)
_download_file(
    "chipseq_samples/ip_2/ip_2.tiny_R1.fastq.gz",
    "workflows/chipseq/data/example_data/chipseq_ip2.fq.gz",
    args.verbose,
)
_download_file(
    "chipseq_samples/ip_3/ip_3.tiny_R1.fastq.gz",
    "workflows/chipseq/data/example_data/chipseq_ip3.fq.gz",
    args.verbose,
)
_download_file(
    "chipseq_samples/ip_4/ip_4.tiny_R1.fastq.gz",
    "workflows/chipseq/data/example_data/chipseq_ip4.fq.gz",
    args.verbose,
)
_download_file(
    "chipseq_samples/input_3/input_3.tiny_R1.fastq.gz",
    "workflows/chipseq/data/example_data/chipseq_input3.fq.gz",
    args.verbose,
)
