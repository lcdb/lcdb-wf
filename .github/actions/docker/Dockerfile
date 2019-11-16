FROM continuumio/miniconda3

COPY requirements.txt requirements.txt
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda create -n lcdb-wf-test --file requirements.txt
