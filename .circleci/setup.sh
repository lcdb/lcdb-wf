#!/bin/bash
set -e

MINICONDA_VER=latest

# Set path
echo "export PATH=$CIRCLE_WORKING_DIRECTORY/miniconda/bin:$PATH" >> $BASH_ENV
source $BASH_ENV

if ! type conda > /dev/null; then
    echo "Setting up conda..."

    # step 1: download and install miniconda
    if [[ $OSTYPE == darwin* ]]; then
        tag="MacOSX"
    elif [[ $OSTYPE == linux* ]]; then
        tag="Linux"
        sudo apt-get install -y fonts-dejavu-core fonts-dejavu-extra freetype2-demos fontconfig
    else
        echo "Unsupported OS: $OSTYPE"
        exit 1
    fi
    curl -L -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-$MINICONDA_VER-$tag-x86_64.sh
    bash miniconda.sh -f -b -p $CIRCLE_WORKING_DIRECTORY/miniconda

    conda config --system --add channels defaults
    conda config --system --add channels bioconda
    conda config --system --add channels conda-forge
    conda config --system --add channels lcdb

    # After SSHing in, for some reason this seems to fix it...
    # conda install -y r-base=3.4.1 bioconductor-genomeinfodbdata bioconductor-annotationhub
    conda create -n lcdb-wf-test -y --file requirements.txt
fi
