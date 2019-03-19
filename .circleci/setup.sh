#!/bin/bash
set -e

WORKSPACE=`pwd`
MINICONDA_VER=latest

# Set path
echo "export PATH=$WORKSPACE/miniconda/bin:$PATH" >> ~/.bashrc
source ~/.bashrc

if ! type conda > /dev/null; then
    echo "Setting up conda..."

    # setup conda if not loaded from cache
    mkdir -p $WORKSPACE

    # step 1: download and install miniconda
    if [[ $OSTYPE == darwin* ]]; then
        tag="MacOSX"
    elif [[ $OSTYPE == linux* ]]; then
        tag="Linux"
    else
        echo "Unsupported OS: $OSTYPE"
        exit 1
    fi
    curl -L -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-$MINICONDA_VER-$tag-x86_64.sh
    bash miniconda.sh -b -p $WORKSPACE/miniconda

    conda config --system --add channels defaults
    conda config --system --add channels bioconda
    conda config --system --add channels conda-forge

    # After SSHing in, for some reason this seems to fix it...
    conda update -y conda
    conda create -n lcdb-wf-test -y --file requirements.txt

fi

