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

cat > ~/.condarc <<EOF
channels:
  - bioconda
  - conda-forge
  - defaults
  - lcdb
default_channels:
  - https://repo.continuum.io/pkgs/free
EOF

# See https://github.com/conda/conda/issues/5536
conda install -y "conda<4.3"

conda --version

echo "Building environment $ENVNAME"
conda create -n $ENVNAME -y --file requirements.txt \
    | grep -v " Time: "

# try pre-caching rnaseq
conda env create -n tmp --file config/envs/R_rnaseq.yaml

source activate $ENVNAME
python ci/get-data.py
