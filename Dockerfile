FROM continuumio/miniconda3

COPY requirements.txt requirements.txt
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
#RUN conda install --file requirements.txt
RUN conda install star
