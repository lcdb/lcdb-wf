#!/bin/bash
set -e

WORKSPACE=`pwd`
MINICONDA_VER=4.3.21

# Set path
echo "export PATH=$WORKSPACE/miniconda/bin:$PATH" >> $BASH_ENV
source $BASH_ENV

cat > ~/.condarc <<EOF
channels:
  - bioconda
  - conda-forge
  - defaults
  - lcdb
default_channels:
  - https://repo.continuum.io/pkgs/free
EOF

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
    conda config --system --add channels conda-forge
    conda config --system --add channels bioconda
    conda config --system --add channels lcdb

    # After SSHing in, for some reason this seems to fix it...
    conda install -y r-base=3.4.1 bioconductor-genomeinfodbdata
    conda create -n lcdb-wf-test -y --file requirements.txt
    conda remove -y r-base
fi

