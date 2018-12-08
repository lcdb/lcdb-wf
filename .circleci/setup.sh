#!/bin/bash
set -e

MINICONDA_VER=latest

# Set path
echo "export PATH=/miniconda/bin:$PATH" >> $BASH_ENV
source $BASH_ENV

cat > ~/.condarc <<EOF
channels:
  - conda-forge
  - bioconda
  - defaults
  - lcdb
default_channels:
  - https://repo.continuum.io/pkgs/main
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
        sudo apt-get install -y fonts-dejavu-core fonts-dejavu-extra freetype2-demos fontconfig
    else
        echo "Unsupported OS: $OSTYPE"
        exit 1
    fi
    curl -L -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-$MINICONDA_VER-$tag-x86_64.sh
    bash miniconda.sh -b -p /miniconda

    #conda config --system --add channels defaults
    #conda config --system --add channels bioconda
    #conda config --system --add channels conda-forge
    #conda config --system --add channels lcdb

    # After SSHing in, for some reason this seems to fix it...
    # conda install -y r-base=3.4.1 bioconductor-genomeinfodbdata bioconductor-annotationhub
    conda create -n lcdb-wf-test -y --file requirements.txt -c conda-forge -c bioconda
fi
