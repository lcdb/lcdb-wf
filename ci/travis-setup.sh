#!/bin/bash
set -eo pipefail
set -x

CONDA_DIR=$HOME/anaconda
ENVNAME=lcdb-wf-test

# Install miniconda if it doesn't already exist.
if [[ ! -e $CONDA_DIR ]]; then
    curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -f -b -p $CONDA_DIR
else
    echo "$CONDA_DIR" exists, so skipping installation of miniconda
fi

export PATH="$CONDA_DIR/bin:$PATH"

# Add channels in the specified order.
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda config --add channels lcdb

# Environment exists; update with requirements. This can happen on a buildkite
# agent.
if conda env list | grep -q $ENVNAME; then
    echo "Environment $ENVNAME exists; updating it"
    conda install -y --file requirements.txt python=3.5
else
    # Otherwise create a new environment (this will happen every time on
    # travis-ci). Don't show the progress for package downloads in the main
    # log, but provide setup.log as a build artifact for buildkite.
    echo "Building environment $ENVNAME"
    conda create -n $ENVNAME -y python=3.5 --file requirements.txt \
        | tee setup.log \
        | grep -v " Time: "
fi

source activate $ENVNAME
python ci/get-data.py

# buildkite expects this, so create if it doesn't exist.
touch setup.log
