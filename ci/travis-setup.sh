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

# Add channels in the specified order.
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda config --add channels lcdb

echo "Building environment $ENVNAME"
conda create -n $ENVNAME -y python=3.5 --file requirements.txt \
    | grep -v " Time: "

source activate $ENVNAME
python ci/get-data.py
