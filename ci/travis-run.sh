#!/bin/bash

set -eo pipefail

if [[ $TRAVIS_BRANCH != "master" && $TRAVIS_PULL_REQUEST == "false" && $TRAVIS_REPO_SLUG == "lcdb/lcdb-wf" ]]
then
    echo ""
    echo "Tests are skipped for non-master-branch pushes to the main repo."
    echo "If you have opened a pull request, please see the full tests for that PR."
    echo ""
    exit 0
fi

case $TYPE in
  references.snakefile)
    cd workflows/references \
    && source activate lcdb-wf-test \
    && snakemake \
    --configfile config/config.yaml \
    --use-conda \
    -j2 -T -k -p -r
    ;;
  rnaseq.snakefile)
    cd workflows/rnaseq \
    && source activate lcdb-wf-test \
    && snakemake \
    --configfile config/config.yaml \
    --use-conda \
    -j2 -T -k -p -r
    ;;
  chipseq.snakefile)
    cd workflows/chipseq \
    && source activate lcdb-wf-test && \
    snakemake \
    --configfile config/config.yaml \
    --use-conda \
    -j2 -T -k -p -r
    ;;
  #pytest)
  #  source activate lcdb-wf-test && py.test wrappers/test -v
  #  ;;
  docs)
    ci/build-docs.sh
    ;;
  figures.snakefile)
    cd workflows/figures \
    && source activate lcdb-wf-test && \
    snakemake \
    --use-conda \
    -j2 -T -k -p -r
    ;;
esac
