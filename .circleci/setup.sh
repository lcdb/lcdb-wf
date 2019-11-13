#!/bin/bash
set -e
if ! conda env list | grep -q "lcdb-wf-test"; then
    echo "Setting up conda..."
    conda config --system --add channels defaults
    conda config --system --add channels bioconda
    conda config --system --add channels conda-forge
    conda create -n lcdb-wf-test -y --file requirements.txt
fi

