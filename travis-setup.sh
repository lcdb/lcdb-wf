#!/bin/bash
set -euo pipefail
set -x

# Sets up travis-ci environment for testing bioconda-utils.
curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/anaconda
export PATH=$HOME/anaconda/bin:$PATH

# Add channels in the specified order.
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda config --add channels lcdb

conda install -y python=3.5
conda install -y --file requirements.txt

python get-data.py
