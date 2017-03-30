#!/bin/bash
set -eo pipefail
set -x

CONDA_DIR=/tmp/ci/ana
ENVNAME=lcdb-wf-test

# buildkite expects this as a build artifact; make sure it exists even if other
# commands fail
touch setup.log

curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -f -b -p $CONDA_DIR

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

    # echo "Environment $ENVNAME exists; updating it"
    # conda install -y --file requirements.txt python=3.5

    # Remove the env. Seems like we're getting strange mkl symlink errors
    conda env remove -y -n $ENVNAME
fi

# Otherwise create a new environment (this will happen every time on
# travis-ci). Don't show the progress for package downloads in the main
# log, but provide setup.log as a build artifact for buildkite.
echo "Building environment $ENVNAME"
conda create -n $ENVNAME -y python=3.5 --file requirements.txt \
    | tee setup.log \
    | grep -v " Time: "

source activate $ENVNAME
python ci/get-data.py
