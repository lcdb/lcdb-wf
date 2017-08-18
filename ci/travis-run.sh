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
    source activate lcdb-wf-test \
      && snakemake -s references.snakefile \
      --configfile config/test_config.yaml \
      --use-conda \
      -j2 -T -k -p -r
    ;;
  rnaseq.snakefile)
    source activate lcdb-wf-test && \
      conda --version && \
      snakemake -s rnaseq.snakefile \
      --configfile config/test_config.yaml \
      --use-conda \
      -j2 -T -k -p -r
    ;;
  chipseq.snakefile)
    source activate lcdb-wf-test && \
      conda --version && \
      snakemake -s chipseq.snakefile \
      --configfile config/test_chipseq_config.yaml \
      --use-conda \
      -j2 -T -k -p -r
    ;;
  pytest)
    source activate lcdb-wf-test && py.test wrappers/test -v
    ;;
  docs)
    ci/build-docs.sh
    ;;
esac
