#!/bin/bash
set -eo pipefail
set -x

if [[ $TRAVIS_BRANCH != "master" && $TRAVIS_PULL_REQUEST == "false" && $TRAVIS_REPO_SLUG == "lcdb/lcdb-wf" ]]
then
    echo ""
    echo "Tests are skipped for non-master-branch pushes to the main repo."
    echo "If you have opened a pull request, please see the full tests for that PR."
    echo ""
    exit 0
fi

CONDA_DIR=/tmp/ci/ana
ENVNAME=lcdb-wf-test

curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -f -b -p $CONDA_DIR

export PATH="$CONDA_DIR/bin:$PATH"

# See https://github.com/conda/conda/issues/5536
conda install -y "conda<4.3"

# Add channels in the specified order.
conda config --add channels conda-forge
conda config --add channels defaults

# Recently bioconda helped migrate a ton of R packages from the `r` channel to
# the `conda-forge` channel. See https://github.com/conda/conda/issues/5536 for
# why the r channel needs to be removed...

# conda config --add channels r
conda config --add channels bioconda
conda config --add channels lcdb

echo "Building environment $ENVNAME"
conda create -n $ENVNAME -y --file requirements.txt \
    | grep -v " Time: "

source activate $ENVNAME
python ci/get-data.py
