#!/bin/bash
set -e

apt install -y locales-all locale
LC_ALL=en_US.utf8
LANG=en_US.utf8
if ! conda env list | grep -q "lcdb-wf-test"; then
    echo "Setting up conda..."
    conda config --system --add channels defaults
    conda config --system --add channels bioconda
    conda config --system --add channels conda-forge
    conda create -n lcdb-wf-test -y --file requirements.txt
fi

