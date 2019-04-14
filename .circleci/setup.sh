#!/bin/bash
set -e

MINICONDA_VER=latest
tag="Linux"

if ! [ -x "$(command -v conda)" ]; then
  apt-get update
  apt-get install -y curl
  curl -L -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-$MINICONDA_VER-$tag-x86_64.sh
  bash miniconda.sh -b -p $CI_PROJECT_DIR/miniconda
  conda update -y conda
fi

export PATH=$CI_PROJECT_DIR/miniconda/bin:$PATH

conda config --system --add channels defaults
conda config --system --add channels bioconda
conda config --system --add channels conda-forge

conda create -n lcdb-wf-test -y --file requirements.txt
