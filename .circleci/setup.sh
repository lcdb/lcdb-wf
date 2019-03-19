#!/bin/bash
set -e

WORKSPACE=`pwd`
MINICONDA_VER=latest
tag="Linux"

apt-get update
apt-get install -y curl

# Set path
curl -L -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-$MINICONDA_VER-$tag-x86_64.sh
bash miniconda.sh -b -p $WORKSPACE/miniconda

set -x
echo "export PATH=$WORKSPACE/miniconda/bin:$PATH" >> $WORKSPACE/.bashrc
source $WORKSPACE/.bashrc
echo $PATH
set +x

conda config --system --add channels defaults
conda config --system --add channels bioconda
conda config --system --add channels conda-forge

# After SSHing in, for some reason this seems to fix it...
conda update -y conda
conda create -n lcdb-wf-test -y --file requirements.txt
