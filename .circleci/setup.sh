#!/bin/bash
set -eu

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

    # step 2: setup channels
    conda config --system --add channels defaults
    conda config --system --add channels conda-forge
    conda config --system --add channels bioconda
    conda config --system --add channels lcdb


    # step 3: install bioconda-utils
    conda install -y --file requirements.txt

    # step 5: cleanup
    conda clean -y --all
fi

# Fetch the master branch for comparison (this can fail locally, if git remote 
# is configured via ssh and this is executed in a container).
git fetch origin +master:master || true
