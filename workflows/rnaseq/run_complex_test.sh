#!/bin/bash
snakemake \
    --configfile ../../test/test_configs/complex-dataset-rnaseq-config.yaml \
    --config sampletable=../../test/test_configs/complex-dataset-rnaseq-sampletable.tsv merged_bigwigs="{}" \
    "$@"
