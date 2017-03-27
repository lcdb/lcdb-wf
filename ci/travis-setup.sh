#!/bin/bash
set -eo pipefail
set -x

# Sets up travis-ci environment for testing bioconda-utils.
curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
mkdir -p /opt/ci/miniconda
bash Miniconda3-latest-Linux-x86_64.sh -f -b -p /opt/ci/miniconda
export PATH=/opt/ci/miniconda/bin:$PATH

# Add channels in the specified order.
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda config --add channels lcdb

ENVNAME=lcdb-wf-test
conda install -y "conda<4.3"
conda env list | grep -q $ENVNAME && conda env remove -y -n $ENVNAME
conda create -n $ENVNAME -y python=3.5 --file requirements.txt | tee setup.log | grep -v " Time: "
source activate $ENVNAME

python ci/get-data.py
