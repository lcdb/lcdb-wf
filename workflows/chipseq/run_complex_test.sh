#!/bin/bash
snakemake \
    --configfile ../../test/test_configs/complex-dataset-chipseq-config.yaml \
    --config \
      sampletable=../../test/test_configs/complex-dataset-chipseq-sampletable.tsv \
    "$@"
